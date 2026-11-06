"""Scalar NumPy example for lesson 06. Run directly to see a worked case."""

from typing import NamedTuple

import numpy as np
from scipy.special import log_ndtr


class ScalarGreeks(NamedTuple):
    vega: float
    theta: float
    rho: float


def _weighted_exp(weight, log_value):
    if weight == 0:
        return 0.0
    return np.copysign(np.exp(np.log(abs(weight)) + log_value), weight)


def greeks_scalar(S, K, r, q, sigma, T, *, kind="call") -> ScalarGreeks:
    """European sensitivities per unit volatility/rate and per calendar year."""
    args = np.asarray([S, K, r, q, sigma, T], dtype=np.float64)
    if args.shape != (6,) or not np.all(np.isfinite(args)):
        raise ValueError("Six finite scalar inputs are required")
    S, K, r, q, sigma, T = args
    if S <= 0 or K <= 0 or sigma < 0 or T < 0:
        raise ValueError("Require S,K>0 and sigma,T>=0")
    if kind not in ("call", "put"):
        raise ValueError("kind must be 'call' or 'put'")
    eta = 1 if kind == "call" else -1

    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        if T == 0:
            if S == K and sigma > 0:
                return ScalarGreeks(0.0, -np.inf, 0.0)
            if eta * (S - K) < 0:
                return ScalarGreeks(0.0, 0.0, 0.0)
            slope = eta * ((q-r)*S if S == K else q*S-r*K)
            theta = min(slope, 0.0) if S == K else slope
            result = ScalarGreeks(0.0, theta, 0.0)
        else:
            # A = S exp(-qT), B = K exp(-rT). Keep their logs until
            # the normal probability/density has been combined with them.
            log_A, log_B = np.log(S) - q*T, np.log(K) - r*T
            relative = (S - K) / K
            log_sk = np.log1p(relative) if abs(relative) < 0.5 else np.log(S) - np.log(K)
            rate_gap = r - q
            carry = rate_gap*T if np.isfinite(rate_gap) else r*T - q*T
            m = log_sk + carry
            if not all(np.isfinite(x) for x in (m, log_A, log_B)):
                raise OverflowError("Black-Scholes log forward or discounts are not representable")
            root_T = np.sqrt(T)
            stddev = sigma * root_T
            if stddev > 0:
                center = m / stddev
            else:
                center = 0.0 if m == 0 else np.copysign(np.inf, m)
            half = 0.5*stddev if np.isfinite(stddev) else (0.5*sigma)*root_T
            d1, d2 = center + half, center - half
            log_Aphi = log_A - 0.5*np.log(2*np.pi) - (0.5*d1)*d1
            log_AN = log_A + log_ndtr(eta*d1)
            log_BN = log_B + log_ndtr(eta*d2)
            vega = np.exp(log_Aphi + 0.5*np.log(T))
            rho = eta * np.exp(log_BN + np.log(T))
            if sigma == 0 and m == 0:
                theta = max(_weighted_exp(-eta*rate_gap, log_A), 0.0)
            else:
                if sigma > 0 and stddev*(1 + abs(center)) <= 1e-3:
                    # Integrate the density over [d2, d1] using its
                    # symmetric expansion around the midpoint.
                    correction = ((center*stddev)**2 - stddev**2)/24
                    log_mass = (np.log(sigma) + 0.5*np.log(T) - (0.5*center)*center
                                - 0.5*np.log(2*np.pi) + np.log1p(correction))
                    carry = eta*_weighted_exp(q*np.expm1(m) - rate_gap, log_BN)
                    carry += _weighted_exp(q, log_A + log_mass)
                else:
                    carry = eta * (_weighted_exp(q, log_AN) - _weighted_exp(r, log_BN))
                theta = carry - _weighted_exp(sigma, log_Aphi - np.log(2*root_T))
            result = ScalarGreeks(vega, theta, rho)

    if not np.all(np.isfinite(result)):
        raise OverflowError("Black-Scholes Greeks are not representable as float64")
    return result


if __name__ == "__main__":
    for kind in ("call", "put"):
        vega, theta, rho = greeks_scalar(100, 100, 0.03, 0.01, 0.20, 1, kind=kind)
        print(f"{kind}: Vega/vol point={vega*0.01:.8f}, "
              f"Theta/year={theta:.8f}, Rho/rate point={rho*0.01:.8f}")
