"""Vectorized European Vega, calendar Theta, and Rho in float64."""

from dataclasses import dataclass

import numpy as np
from scipy.special import log_ndtr

from black_scholes import LOG_SQRT_2PI
from bs_batch import _validated_inputs


@dataclass(frozen=True)
class GreeksResult:
    vega: np.ndarray
    theta: np.ndarray
    rho: np.ndarray


def _weighted_exp(coefficient, log_value):
    """Multiply before exponentiating so a tiny tail can survive rescaling."""
    return np.where(
        coefficient == 0, 0.0,
        np.sign(coefficient) * np.exp(np.log(np.abs(coefficient)) + log_value),
    )


def evaluate_greeks(S, K, r, q, sigma, T, *, kind="call") -> GreeksResult:
    """Broadcast inputs in the same order as evaluate_call.

    Vega and Rho are per unit change in sigma and r. Theta is per year.
    At positive-volatility ATM expiry, Theta is the extended limit -inf.
    At zero-volatility forward ATM, Vega is a right derivative and Rho
    uses the midpoint of the left and right derivatives. See docs/bs_greeks.md.
    """
    if kind not in ("call", "put"):
        raise ValueError("kind must be 'call' or 'put'")
    sign = 1 if kind == "call" else -1
    values = _validated_inputs(S, K, r, q, sigma, T)
    shape = values[0].shape
    s, k, rate, yield_, vol, t = (value.ravel() for value in values)
    vega, theta, rho = (np.zeros_like(s) for _ in range(3))
    live = t > 0
    expiry_atm = ~live & (s == k)
    singular = expiry_atm & (vol > 0)
    with np.errstate(over="ignore", under="ignore", divide="ignore", invalid="ignore"):
        expired_itm = ~live & (sign * (s - k) > 0)
        edge = expired_itm | (expiry_atm & ~singular)
        expiry_carry = sign * (yield_[edge] * s[edge] - rate[edge] * k[edge])
        atm_edge = expiry_atm & edge
        expiry_carry[expiry_atm[edge]] = sign * (yield_[atm_edge] - rate[atm_edge]) * s[atm_edge]
        theta[edge] = np.where(expired_itm[edge], expiry_carry, np.minimum(expiry_carry, 0))
        theta[singular] = -np.inf

        s, k, r, q, vol, t = (x[live] for x in (s, k, rate, yield_, vol, t))
        root_t = np.sqrt(t)
        log_s, log_k = np.log(s), np.log(k)
        relative = (s - k) / k
        log_sk = np.where(np.abs(relative) < 0.5, np.log1p(relative), log_s - log_k)
        carry_rate = r - q
        rt, qt = r * t, q * t
        m = log_sk + np.where(np.isfinite(carry_rate), carry_rate * t, rt - qt)
        log_a, log_b = log_s - qt, log_k - rt
        if any(not np.all(np.isfinite(x)) for x in (m, log_a, log_b)):
            raise OverflowError("Black-Scholes log forward or discounts are not representable")

        stddev = vol * root_t
        # Divide the log-moneyness first. Form d2 directly instead of
        # subtracting two large numbers when standard deviation is large.
        limit = np.where(m == 0, 0.0, np.copysign(np.inf, m))
        scaled_m = np.divide(m, stddev, out=limit, where=stddev > 0)
        half = np.where(np.isfinite(stddev), 0.5 * stddev, (0.5 * vol) * root_t)
        d1, d2 = scaled_m + half, scaled_m - half
        log_density = log_a - LOG_SQRT_2PI - (0.5 * d1) * d1
        log_spot_leg = log_a + log_ndtr(sign * d1)
        log_strike_leg = log_b + log_ndtr(sign * d2)
        carry = sign * (
            _weighted_exp(q, log_spot_leg) - _weighted_exp(r, log_strike_leg)
        )
        # Evaluate the probability of a narrow interval directly, so the
        # carry survives cancellation between the two discounted legs.
        narrow = (vol > 0) & (stddev * (1 + np.abs(scaled_m)) <= 1e-3)
        c, tau = scaled_m[narrow], stddev[narrow]
        correction = ((c*tau)**2 - tau**2) / 24
        log_mass = (np.log(vol[narrow]) + 0.5*np.log(t[narrow])
                    - (0.5*c)*c - LOG_SQRT_2PI + np.log1p(correction))
        carry[narrow] = (
            sign * _weighted_exp(q[narrow]*np.expm1(m[narrow]) - carry_rate[narrow], log_strike_leg[narrow])
            + _weighted_exp(q[narrow], log_a[narrow] + log_mass)
        )
        vega[live] = np.exp(log_density + 0.5 * np.log(t))
        rho[live] = sign * np.exp(log_strike_leg + np.log(t))
        diffusion = _weighted_exp(vol, log_density - np.log(2 * root_t))
        # At the deterministic kink the calendar-right derivative is the
        # positive part of the slope of signed discounted intrinsic value.
        kink = (vol == 0) & (m == 0)
        kink_carry = _weighted_exp(-sign * carry_rate[kink], log_a[kink])
        live_theta = carry - diffusion
        live_theta[kink] = np.maximum(kink_carry, 0.0)
        theta[live] = live_theta

    if not (np.all(np.isfinite(vega)) and np.all(np.isfinite(rho))
            and np.all(np.isfinite(theta[~singular]))):
        raise OverflowError("Black-Scholes Greeks are not representable as float64")
    return GreeksResult(vega.reshape(shape), theta.reshape(shape), rho.reshape(shape))
