"""Numerical and API contracts shared by the scalar and vector Greeks."""

import csv
import math
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python_prototype"))
from bs_batch import evaluate_call
from bs_greeks import evaluate_greeks
from greeks_prototype import greeks_scalar

INPUTS = ("S", "K", "r", "q", "sigma", "T")
FIELDS = ("vega", "theta", "rho")


class GreeksTests(unittest.TestCase):
    def test_shared_high_precision_fixture(self):
        with (ROOT / "testdata" / "bs_greeks_golden.csv").open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        for kind in ("call", "put"):
            selected = [row for row in rows if row["kind"] == kind]
            self.assertGreaterEqual(len(selected), 20)
            args = np.array([[float(row[key]) for key in INPUTS] for row in selected])
            batch = evaluate_greeks(*args.T, kind=kind)
            for i, row in enumerate(selected):
                with self.subTest(case=row["case"], kind=kind):
                    expected = [float(row[key]) for key in FIELDS]
                    scalar = greeks_scalar(*args[i], kind=kind)
                    vector = [getattr(batch, key)[i] for key in FIELDS]
                    np.testing.assert_allclose(scalar, expected, rtol=2e-10, atol=0)
                    np.testing.assert_allclose(vector, expected, rtol=2e-10, atol=0)

    def test_broadcast_empty_and_unchanged_input(self):
        spots = np.array([[90.0], [110.0]])
        strikes = np.array([80.0, 100.0, 120.0])
        for kind in ("call", "put"):
            grid = evaluate_greeks(spots, strikes, -.01, .02, .2, 1, kind=kind)
            one = evaluate_greeks(100, 100, .01, .02, .2, 1, kind=kind)
            empty = evaluate_greeks(np.empty((0, 1)), strikes, 0, 0, .2, 1, kind=kind)
            for key in FIELDS:
                self.assertEqual(getattr(one, key).shape, ())
                self.assertEqual(getattr(empty, key).shape, (0, 3))
                self.assertEqual(getattr(grid, key).shape, (2, 3))
            for i, s in enumerate(spots[:, 0]):
                for j, k in enumerate(strikes):
                    expected = greeks_scalar(s, k, -.01, .02, .2, 1, kind=kind)
                    np.testing.assert_allclose([getattr(grid, key)[i, j] for key in FIELDS], expected)
        np.testing.assert_array_equal(spots, [[90], [110]])

    def test_invalid_and_overflow(self):
        for fn in (evaluate_greeks, greeks_scalar):
            for args in ((0, 100, 0, 0, .2, 1), (100, 0, 0, 0, .2, 1),
                         (100, 100, math.inf, 0, .2, 1), (100, 100, 0, math.nan, .2, 1),
                         (100, 100, 0, 0, -.2, 1), (100, 100, 0, 0, .2, -1)):
                with self.subTest(fn=fn.__name__, args=args), self.assertRaises(ValueError):
                    fn(*args)
            with self.assertRaises(ValueError):
                fn(100, 100, 0, 0, .2, 1, kind="straddle")
            with self.assertRaises(OverflowError):
                fn(100, 100, -1000, -1000, .2, 1)
        with self.assertRaises(ValueError):
            evaluate_greeks([1, 2], [1, 2, 3], 0, 0, .2, 1)

    def test_calendar_theta_and_rate_hold_fixed(self):
        args = [108., 101., -.017, .024, .27, 1.3]

        def price(xs, kind):
            s, k, r, q, _, t = xs
            call = float(evaluate_call(*xs).price)
            return call if kind == "call" else call - s*math.exp(-q*t) + k*math.exp(-r*t)

        for kind in ("call", "put"):
            analytical = evaluate_greeks(*args, kind=kind)
            for field, index, sign in (("vega", 4, 1), ("theta", 5, -1), ("rho", 2, 1)):
                with self.subTest(field=field, kind=kind):
                    errors = []
                    for h in (.002, .001, 1e-5):
                        hi, lo = args.copy(), args.copy()
                        hi[index] += h
                        lo[index] -= h
                        fd = sign * (price(hi, kind)-price(lo, kind)) / (2*h)
                        errors.append(abs(fd - float(getattr(analytical, field))))
                    self.assertGreater(errors[0], 3.8 * errors[1])
                    self.assertLess(errors[-1], 1e-7)

    def test_pde_gamma_and_put_call_identities(self):
        s, k, r, q, vol, t = 115., 100., -.025, .055, .32, 1.7
        a, b = s*math.exp(-q*t), k*math.exp(-r*t)
        base = evaluate_call(s, k, r, q, vol, t)
        call = evaluate_greeks(s, k, r, q, vol, t)
        put = evaluate_greeks(s, k, r, q, vol, t, kind="put")
        np.testing.assert_allclose(call.vega, s*s*vol*t*base.gamma, rtol=2e-14)
        np.testing.assert_equal(call.vega, put.vega)
        np.testing.assert_allclose(call.rho-put.rho, t*b, rtol=2e-14)
        np.testing.assert_allclose(call.theta-put.theta, q*a-r*b, rtol=2e-14)
        for greeks, price, delta in ((call, base.price, base.delta),
                                     (put, base.price-a+b, base.delta-math.exp(-q*t))):
            residual = greeks.theta + .5*vol**2*s**2*base.gamma + (r-q)*s*delta - r*price
            self.assertAlmostEqual(float(residual), 0.0, places=12)

    def test_deterministic_kink_one_sided_time(self):
        # S/K=1/2 and r-q=log(2) give exactly zero computed log-moneyness.
        rate = math.log(2)
        for fn in (evaluate_greeks, greeks_scalar):
            call = fn(1, 2, rate, 0, 0, 1, kind="call")
            put = fn(1, 2, rate, 0, 0, 1, kind="put")
            self.assertAlmostEqual(float(call.vega), 1/math.sqrt(2*math.pi), places=14)
            self.assertAlmostEqual(float(call.rho), .5, places=14)
            self.assertAlmostEqual(float(put.rho), -.5, places=14)
            self.assertEqual(float(call.theta), 0.0)
            self.assertAlmostEqual(float(put.theta), rate, places=14)


if __name__ == "__main__":
    unittest.main()
