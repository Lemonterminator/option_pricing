#include "bs_batch.h"

#include "black_scholes.h"
#include "complex_step.h"

#include <algorithm>
#include <cmath>
#include <execution>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace black_scholes {
namespace {

bool valid_inputs(double S, double K, double r, double q, double sigma, double T) noexcept {
    return std::isfinite(S) && std::isfinite(K) && std::isfinite(r) && std::isfinite(q)
        && std::isfinite(sigma) && std::isfinite(T) && S > 0 && K > 0 && sigma >= 0 && T >= 0;
}

void validate(const CallInputs& input) {
    const auto n = input.S.size();
    if (input.K.size() != n || input.r.size() != n || input.q.size() != n ||
        input.sigma.size() != n || input.T.size() != n) {
        throw std::invalid_argument("All batch inputs must have the same length");
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (!valid_inputs(input.S[i], input.K[i], input.r[i], input.q[i], input.sigma[i], input.T[i])) {
            throw std::invalid_argument("Require finite S,K>0, finite r,q, and finite sigma,T>=0");
        }
    }
}

constexpr double log_sqrt_2pi = 0.91893853320467274178;
constexpr double inf = std::numeric_limits<double>::infinity();

double option_sign(OptionKind kind) {
    switch (kind) {
        case OptionKind::Call: return 1.0;
        case OptionKind::Put: return -1.0;
    }
    throw std::invalid_argument("Unknown option kind");
}

double log_cdf(double z) noexcept {
    if (z >= -35) return std::log(Phi_real(z));
    // Normal tail via the Mills expansion. At |z| >= 35 the first
    // omitted term after eight corrections is below double precision.
    const double inv_sq = (1 / z) / z;
    double term = 1;
    double correction = 0;
    for (int n = 1; n <= 8; ++n) {
        term *= -(2 * n - 1) * inv_sq;
        correction += term;
    }
    return -(0.5 * z) * z - log_sqrt_2pi - std::log(-z) + std::log1p(correction);
}

double weighted_exp(double weight, double log_value) noexcept {
    if (weight == 0) return 0;
    return std::copysign(std::exp(std::log(std::abs(weight)) + log_value), weight);
}

ScalarGreeks greeks_unchecked(double S, double K, double r, double q, double sigma, double T,
    double sign) noexcept {
    if (T == 0) {
        if (S == K && sigma > 0) return {0, -inf, 0};
        if (sign * (S - K) < 0) return {0, 0, 0};
        const double slope = sign * (S == K ? (q-r)*S : q*S - r*K);
        return {0, S == K ? std::min(slope, 0.0) : slope, 0};
    }

    const double log_s = std::log(S);
    const double log_k = std::log(K);
    const double relative = (S - K) / K;
    const double log_sk = std::abs(relative) < 0.5 ? std::log1p(relative) : log_s - log_k;
    const double rate_gap = r - q;
    const double rt = r * T;
    const double qt = q * T;
    const double m = log_sk + (std::isfinite(rate_gap) ? rate_gap*T : rt - qt);
    const double log_a = log_s - qt;
    const double log_b = log_k - rt;
    if (!(std::isfinite(m) && std::isfinite(log_a) && std::isfinite(log_b))) {
        return {inf, inf, inf};  // Checked outside the execution-policy callback.
    }

    const double root_t = std::sqrt(T);
    const double stddev = sigma * root_t;
    const double center = stddev > 0 ? m/stddev : (m == 0 ? 0.0 : std::copysign(inf, m));
    const double half = std::isfinite(stddev) ? 0.5*stddev : (0.5*sigma)*root_t;
    const double d1 = center + half;
    const double d2 = center - half;
    const double log_density = log_a - log_sqrt_2pi - (0.5*d1)*d1;
    const double log_spot_leg = log_a + log_cdf(sign*d1);
    const double log_strike_leg = log_b + log_cdf(sign*d2);
    const double vega = std::exp(log_density + 0.5*std::log(T));
    const double rho = sign * std::exp(log_strike_leg + std::log(T));
    if (sigma == 0 && m == 0) {
        return {vega, std::max(weighted_exp(-sign*rate_gap, log_a), 0.0), rho};
    }
    double carry;
    if (sigma > 0 && stddev*(1 + std::abs(center)) <= 1e-3) {
        const double c_tau = center*stddev;
        const double correction = (c_tau*c_tau - stddev*stddev) / 24;
        const double log_mass = std::log(sigma) + 0.5*std::log(T) - (0.5*center)*center
            - log_sqrt_2pi + std::log1p(correction);
        carry = sign*weighted_exp(q*std::expm1(m) - rate_gap, log_strike_leg)
            + weighted_exp(q, log_a + log_mass);
    } else {
        carry = sign * (weighted_exp(q, log_spot_leg) - weighted_exp(r, log_strike_leg));
    }
    const double diffusion = weighted_exp(sigma, log_density - std::log(2*root_t));
    return {vega, carry - diffusion, rho};
}

void check_greeks(const ScalarGreeks& result, bool expiry_singularity) {
    if (!std::isfinite(result.vega) || !std::isfinite(result.rho)
        || !(std::isfinite(result.theta) || (expiry_singularity && result.theta == -inf))) {
        throw std::overflow_error("Black-Scholes Greeks are not representable as double");
    }
}

} // namespace

ScalarGreeks evaluate_greeks(double S, double K, double r, double q, double sigma, double T,
    OptionKind kind) {
    if (!valid_inputs(S, K, r, q, sigma, T)) {
        throw std::invalid_argument("Require finite S,K>0, finite r,q, and finite sigma,T>=0");
    }
    const auto result = greeks_unchecked(S, K, r, q, sigma, T, option_sign(kind));
    check_greeks(result, T == 0 && S == K && sigma > 0);
    return result;
}

GreeksResult evaluate_greeks(const CallInputs& input, OptionKind kind, BatchBackend backend) {
    validate(input);
    const double sign = option_sign(kind);
    const auto n = input.S.size();
    GreeksResult result{std::vector<double>(n), std::vector<double>(n), std::vector<double>(n)};
    const auto evaluate_at = [&](std::size_t i) noexcept {
        const auto row = greeks_unchecked(
            input.S[i], input.K[i], input.r[i], input.q[i], input.sigma[i], input.T[i], sign);
        result.vega[i] = row.vega;
        result.theta[i] = row.theta;
        result.rho[i] = row.rho;
    };
    switch (backend) {
        case BatchBackend::Serial:
            for (std::size_t i = 0; i < n; ++i) evaluate_at(i);
            break;
        case BatchBackend::ParUnseq: {
            std::vector<std::size_t> indices(n);
            std::iota(indices.begin(), indices.end(), std::size_t{0});
            std::for_each(std::execution::par_unseq, indices.begin(), indices.end(), evaluate_at);
            break;
        }
        default: throw std::invalid_argument("Unknown batch backend");
    }
    for (std::size_t i = 0; i < n; ++i) {
        check_greeks({result.vega[i], result.theta[i], result.rho[i]},
            input.T[i] == 0 && input.S[i] == input.K[i] && input.sigma[i] > 0);
    }
    return result;
}

CallResult evaluate_call(const CallInputs& input, BatchBackend backend) {
    // Validate before entering par_unseq: an exception from a parallel element
    // callback would invoke std::terminate under standard execution policies.
    validate(input);
    CallResult result;
    const auto n = input.S.size();
    result.price.resize(n);
    result.delta.resize(n);
    result.gamma.resize(n);

    if (backend == BatchBackend::Serial) {
        for (std::size_t i = 0; i < n; ++i) {
            result.price[i] = bs_price_call(input.S[i], input.K[i], input.r[i], input.q[i], input.sigma[i], input.T[i]);
            result.delta[i] = delta_call(input.S[i], input.K[i], input.r[i], input.q[i], input.sigma[i], input.T[i]);
            result.gamma[i] = Gamma(input.S[i], input.K[i], input.r[i], input.q[i], input.sigma[i], input.T[i]);
        }
    } else if (backend == BatchBackend::ParUnseq) {
        bs_price_call_par(input.S, input.K, input.r, input.q, input.sigma, input.T, result.price);
        delta_call_par(input.S, input.K, input.r, input.q, input.sigma, input.T, result.delta);
        Gamma_par(input.S, input.K, input.r, input.q, input.sigma, input.T, result.gamma);
    } else {
        throw std::invalid_argument("Unknown batch backend");
    }

    for (std::size_t i = 0; i < n; ++i) {
        if (input.sigma[i] == 0.0 || input.T[i] == 0.0) {
            const double a = input.S[i] * std::exp(-input.q[i] * input.T[i]);
            const double b = input.K[i] * std::exp(-input.r[i] * input.T[i]);
            result.price[i] = std::max(a - b, 0.0);
            result.delta[i] = std::exp(-input.q[i] * input.T[i]) *
                (a > b ? 1.0 : a == b ? 0.5 : 0.0);
            result.gamma[i] = 0.0;
        }
        if (!std::isfinite(result.price[i]) || !std::isfinite(result.delta[i]) ||
            !std::isfinite(result.gamma[i])) {
            throw std::overflow_error("Black-Scholes result is not representable as double");
        }
    }
    return result;
}

CallResult evaluate_call_complex_step(const CallInputs& input, const std::vector<double>& h, int order) {
    if (order < 2 || order > 6) throw std::invalid_argument("order must be 2, 3, 4, 5, 6");
    if (h.size() != input.S.size()) throw std::invalid_argument("Step array must match batch length");
    for (double step : h) {
        if (!std::isfinite(step) || step <= 0) throw std::invalid_argument("Steps must be positive and finite");
    }
    auto result = evaluate_call(input);
    for (std::size_t i = 0; i < input.S.size(); ++i) {
        if (input.sigma[i] == 0 || input.T[i] == 0) continue;
        result.delta[i] = delta_complex_step(input.S[i], input.K[i], input.r[i], input.q[i], input.sigma[i], input.T[i], h[i], order);
        result.gamma[i] = gamma_45_deg(input.S[i], input.K[i], input.r[i], input.q[i], input.sigma[i], input.T[i], h[i], order);
    }
    return result;
}

CallResult evaluate_call_complex_step(const CallInputs& input, double h, int order) {
    if (!std::isfinite(h) || h <= 0) throw std::invalid_argument("h must be positive and finite");
    return evaluate_call_complex_step(input, std::vector<double>(input.S.size(), h), order);
}

} // namespace black_scholes
