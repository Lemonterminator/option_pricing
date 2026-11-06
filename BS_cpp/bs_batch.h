/** Validated batch facade over the existing Black-Scholes call API.
 *
 * Backends intentionally share an input/output contract. Parallel requests
 * use the current par_unseq functions; the standard does not guarantee either
 * multiple worker threads or SIMD instructions on a particular build.
 */
#pragma once

#include <vector>

namespace black_scholes {

struct CallInputs {
    const std::vector<double>& S;
    const std::vector<double>& K;
    const std::vector<double>& r;
    const std::vector<double>& q;
    const std::vector<double>& sigma;
    const std::vector<double>& T;
};

struct CallResult {
    std::vector<double> price;
    std::vector<double> delta;
    std::vector<double> gamma;
};

enum class BatchBackend { Serial, ParUnseq };

enum class OptionKind { Call, Put };

struct ScalarGreeks {
    double vega;
    double theta;
    double rho;
};

struct GreeksResult {
    std::vector<double> vega;
    std::vector<double> theta;
    std::vector<double> rho;
};

// Vega and Rho per unit change in sigma and r, calendar Theta per year.
// At ATM expiry with positive volatility, Theta is the extended limit -inf.
// See docs/bs_greeks.md for the deterministic kink conventions.
ScalarGreeks evaluate_greeks(double S, double K, double r, double q, double sigma, double T,
    OptionKind kind = OptionKind::Call);
GreeksResult evaluate_greeks(const CallInputs& input, OptionKind kind = OptionKind::Call,
    BatchBackend backend = BatchBackend::Serial);

// Throws invalid_argument for mismatched, nonfinite or out-of-domain inputs;
// throws overflow_error if an otherwise accepted input has an unrepresentable
// result. Degenerate sigma=0 or T=0 uses intrinsic price, gamma=0, and the
// explicit midpoint delta convention at the intrinsic kink.
CallResult evaluate_call(const CallInputs& input, BatchBackend backend = BatchBackend::Serial);

// Real price plus CS Delta/45-degree Gamma. Analytic degenerate conventions
// are retained. h is an absolute spot step; order is the CDF Taylor power.
CallResult evaluate_call_complex_step(const CallInputs& input,
    const std::vector<double>& h, int order = 3);
CallResult evaluate_call_complex_step(const CallInputs& input,
    double h = 1e-5, int order = 3);

} // namespace black_scholes
