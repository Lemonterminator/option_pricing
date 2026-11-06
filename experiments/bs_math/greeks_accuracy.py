"""Reproduce the shared golden CSV and compare Greeks to mpmath and QuantLib."""

import argparse
import csv
import io
import json
import math
import platform
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Python_prototype"))
from bs_greeks import evaluate_greeks
from bs_batch import evaluate_call
from greeks_prototype import greeks_scalar

INPUTS = ("S", "K", "r", "q", "sigma", "T")
GREEKS = ("vega", "theta", "rho")
KINDS = ("call", "put")
DPS = 80
SEED = 20261106
GOLDEN = (
    ("atm", (100, 100, 0, 0, .2, 1)),
    ("carry", (100, 100, .03, .01, .2, 1)),
    ("negative_rate", (90, 105, -.04, .025, .35, 2)),
    ("positive_theta_call", (160, 100, .01, .15, .06, 1)),
    ("positive_theta_put", (60, 100, .08, .01, .08, 1)),
    ("short_atm", (100, 100, .03, .01, .2, 1e-10)),
    ("low_vol_atm", (100, 100, .03, .03, 1e-10, 1)),
    ("large_stddev", (100, 100, .02, .01, 10, 8)),
    ("deep_call_tail", (25, 100, .01, .02, .11, 1)),
    ("deep_put_tail", (400, 100, .01, .02, .11, 1)),
    ("rescaled_call_tail", (1e196, 3e199, .02, .01, .2, 1)),
    ("rescaled_put_tail", (3e199, 1e196, .02, .01, .2, 1)),
    ("tiny_notional", (1e-200, 1.1e-200, -.02, .01, .2, 2)),
    ("deterministic_itm", (120, 100, .03, .01, 0, 1)),
    ("deterministic_otm", (80, 100, -.03, .02, 0, 1)),
    ("deterministic_atm", (100, 100, .03, .03, 0, 1)),
    ("expiry_itm", (120, 100, .03, .01, .2, 0)),
    ("expiry_otm", (80, 100, -.03, .02, .2, 0)),
    ("expiry_atm", (100, 100, .03, .01, .2, 0)),
    ("deterministic_expiry_atm", (100, 100, .03, .01, 0, 0)),
)


def mp_reference(args, kind):
    """Differentiate high-precision prices with all other inputs held fixed."""
    s, k, r, q, vol, t = (mp.mpf(float(x)) for x in args)
    eta = 1 if kind == "call" else -1
    a, b = s*mp.exp(-q*t), k*mp.exp(-r*t)
    carry = eta*(q*a-r*b)
    if t == 0:
        if s == k and vol > 0:
            return 0, -mp.inf, 0
        theta = carry if eta*(s-k) > 0 else min(carry, 0) if s == k else 0
        return 0, theta, 0
    if vol == 0:
        if a == b:
            return a*mp.sqrt(t/(2*mp.pi)), max(carry, 0), eta*t*b/2
        return (0, carry, eta*t*b) if eta*(a-b) > 0 else (0, 0, 0)

    # Differentiate the OTM premium and use parity for the ITM option.
    # This avoids subtracting nearly equal intrinsic values in a tiny Vega.
    use_parity = eta*(a-b) > 0
    price_sign = -eta if use_parity else eta

    def price(r_, v_, t_):
        a_, b_ = s*mp.exp(-q*t_), k*mp.exp(-r_*t_)
        tau = v_*mp.sqrt(t_)
        d1 = (mp.log(s/k)+(r_-q)*t_)/tau + tau/2
        d2 = d1-tau
        n1 = mp.erfc(-price_sign*d1/mp.sqrt(2))/2
        n2 = mp.erfc(-price_sign*d2/mp.sqrt(2))/2
        return price_sign*(a_*n1-b_*n2)

    vega = mp.diff(lambda v_: price(r, v_, t), vol)
    theta = -mp.diff(lambda t_: price(r, vol, t_), t)
    rho = mp.diff(lambda r_: price(r_, vol, t), r)
    return vega, theta + carry if use_parity else theta, rho + eta*t*b if use_parity else rho


def fixtures(cases):
    with mp.workdps(DPS):
        return [(name, kind, *args, *(mp.nstr(x, 45) for x in mp_reference(args, kind)))
                for name, args in cases for kind in KINDS]


def write_fixture(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("case", "kind", *INPUTS, *GREEKS))
        writer.writerows(rows)


def sampled_cases():
    rng = np.random.default_rng(SEED)
    for i in range(384):
        strike = rng.uniform(30, 220)
        yield f"regular_{i:03}", (
            strike*math.exp(rng.uniform(-.6, .6)), strike,
            rng.uniform(-.06, .15), rng.uniform(-.04, .14),
            rng.uniform(.08, .9), rng.uniform(.06, 12),
        )
    for t in (1e-12, 1e-7, 1/365, .25, 5):
        for vol in (1e-6, 1e-3, .2, 2):
            for center in (-40, -12, -2, 0, 2, 12, 40):
                # Target d1/d2's midpoint in standard-deviation units.
                k, r, q = 100, .03, .01
                spot = k*math.exp(center*vol*math.sqrt(t)-(r-q)*t)
                yield f"stress_{t:g}_{vol:g}_{center}", (spot, k, r, q, vol, t)


def quantlib_values(ql, args, kind):
    s, k, r, q, vol, t = args
    payoff = ql.PlainVanillaPayoff(ql.Option.Call if kind == "call" else ql.Option.Put, k)
    calculator = ql.BlackCalculator(payoff, s*math.exp((r-q)*t), vol*math.sqrt(t), math.exp(-r*t))
    return calculator.vega(t), calculator.theta(s, t), calculator.rho(t)


def scale(args, greek):
    s, k, r, q, vol, t = args
    a, b = s*math.exp(-q*t), k*math.exp(-r*t)
    if greek == "vega":
        return a*math.sqrt(t)
    if greek == "rho":
        return b*t
    return abs(q*a) + abs(r*b) + (a*vol/(2*math.sqrt(t)) if t else 0.0)


def finite_difference_sweep(output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    scenarios = {"1 year": (100, 100, .03, .01, .2, 1),
                 "1 day": (100, 100, .03, .01, .2, 1/365)}
    bumps = np.geomspace(.1, 1e-15, 57)
    figure, axes = plt.subplots(1, 3, figsize=(10.2, 3.25), layout="constrained")
    with (output / "finite_difference.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("scenario", "greek", "h_rel", "h", "estimate", "reference", "abs_error", "status"))
        for label, args in scenarios.items():
            with mp.workdps(DPS):
                refs = [float(x) for x in mp_reference(args, "call")]
            for ax, greek, index, sign, ref in zip(axes, GREEKS, (4, 5, 2), (1, -1, 1), refs):
                x = args[index]
                steps = bumps * abs(x)
                errors = []
                for h_rel, h in zip(bumps, steps):
                    hi, lo = list(args), list(args)
                    hi[index], lo[index] = x+h, x-h
                    if hi[index] == x or lo[index] == x:
                        writer.writerow((label, greek, h_rel, h, "", ref, "", "collapsed bump"))
                        errors.append(np.nan)
                        continue
                    estimate = sign * float(evaluate_call(*hi).price - evaluate_call(*lo).price) / (2*h)
                    error = abs(estimate-ref)
                    writer.writerow((label, greek, h_rel, h, estimate, ref, error, "compared"))
                    errors.append(error if error > 0 else np.nan)
                ax.loglog(bumps, errors, marker=".", markersize=3, linewidth=1, label=label)
                ax.set(title=greek.title(), xlabel="Relative bump h / |x|")
                ax.grid(which="major", alpha=.24)
    for ax in axes:
        ax.invert_xaxis()
    axes[0].set_ylabel("Absolute error vs mpmath")
    axes[-1].legend(loc="best", frameon=False)
    figure.savefig(output / "finite_difference.png", dpi=180)
    plt.close(figure)


def compare(executable, output):
    import QuantLib as ql

    rows = fixtures((*GOLDEN, *sampled_cases()))
    output.mkdir(parents=True, exist_ok=True)
    fixture_path = output / "reference.csv"
    write_fixture(fixture_path, rows)
    cpp = subprocess.run([str(executable), str(fixture_path), "--dump"], check=True,
                         capture_output=True, text=True, timeout=90)
    methods = defaultdict(dict)
    for row in csv.DictReader(io.StringIO(cpp.stdout)):
        methods[row["method"]][row["case"], row["kind"]] = tuple(float(row[g]) for g in GREEKS)
    if set(methods) != {"cpp_scalar", "cpp_batch", "cpp_par_unseq"}:
        raise RuntimeError("The C++ probe did not emit all three backends")
    for kind in KINDS:
        selected = [row for row in rows if row[1] == kind]
        batch = evaluate_greeks(*np.asarray([row[2:8] for row in selected]).T, kind=kind)
        for i, (name, _, *fields) in enumerate(selected):
            args = fields[:6]
            key = name, kind
            methods["numpy_scalar"][key] = greeks_scalar(*args, kind=kind)
            methods["numpy_vector"][key] = tuple(getattr(batch, greek)[i] for greek in GREEKS)
            if args[-2] == 0 or args[-1] == 0:
                methods["quantlib"][key] = "boundary convention"
                continue
            try:
                methods["quantlib"][key] = quantlib_values(ql, args, kind)
            except (RuntimeError, OverflowError) as exc:
                methods["quantlib"][key] = str(exc)

    summary = defaultdict(lambda: dict(compared=0, unavailable=0, max_abs=0.0,
                                       max_rel_normal=0.0, max_scaled=0.0,
                                       zeroed_normal=0, worst_rel_case="", worst_scaled_case=""))
    failures = []
    with (output / "comparison.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("case", "kind", "group", "greek", "method", "reference",
                         "actual", "abs_error", "rel_error", "scaled_error", "status"))
        for name, kind, *fields in rows:
            args, refs = fields[:6], fields[6:]
            boundary = args[-2] == 0 or args[-1] == 0
            group = ("boundary" if boundary else "regular" if name.startswith("regular_") else
                     "stress" if name.startswith("stress_") else "golden")
            for method, values in methods.items():
                result = values[name, kind]
                for i, (greek, ref) in enumerate(zip(GREEKS, refs)):
                    record = summary[group, method, greek]
                    ref64 = float(ref)
                    actual = math.nan if isinstance(result, str) else float(result[i])
                    if method == "quantlib" and not math.isfinite(actual):
                        record["unavailable"] += 1
                        if group == "regular":
                            failures.append((name, kind, greek, method, "unavailable"))
                        writer.writerow((name, kind, group, greek, method, ref, "", "", "", "",
                                         result if isinstance(result, str) else "nonfinite QuantLib result"))
                        continue
                    abs_error = 0.0 if actual == ref64 else abs(actual-ref64)
                    unit = scale(args, greek)
                    scaled = (abs_error/unit if unit else abs_error) if math.isfinite(abs_error) else math.inf
                    relative = abs_error/abs(ref64) if math.isfinite(ref64) and ref64 else None
                    record["compared"] += 1
                    record["max_abs"] = max(record["max_abs"], abs_error)
                    if scaled > record["max_scaled"]:
                        record["max_scaled"], record["worst_scaled_case"] = scaled, f"{name}/{kind}"
                    if relative is not None and abs(ref64) >= sys.float_info.min:
                        if relative > record["max_rel_normal"]:
                            record["max_rel_normal"], record["worst_rel_case"] = relative, f"{name}/{kind}"
                        if actual == 0:
                            record["zeroed_normal"] += 1
                    status = "compared"
                    if (method != "quantlib" or group == "regular") and (
                        not math.isfinite(scaled) or scaled > 5e-11
                    ):
                        status = "outside scaled tolerance"
                        failures.append((name, kind, greek, method, scaled))
                    writer.writerow((name, kind, group, greek, method, ref, actual, abs_error,
                                     relative, scaled, status))

    records = [dict(group=g, method=m, greek=k, **v) for (g, m, k), v in sorted(summary.items())]
    payload = {
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "scipy": scipy.__version__, "mpmath": mp.__version__, "quantlib": ql.__version__,
                        "platform": platform.system(), "machine": platform.machine()},
        "reference_dps": DPS, "seed": SEED, "case_count": len(rows), "scaled_tolerance": 5e-11,
        "summary": records, "failures": failures,
    }
    with (output / "accuracy_summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=records[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    (output / "results.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    finite_difference_sweep(output)
    print(f"{len(rows)} option cases; QuantLib {ql.__version__}; {len(failures)} tolerance failures")
    return bool(failures)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-golden", action="store_true")
    parser.add_argument("--cpp", type=Path, help="Path to the built bs_greeks_contract_test")
    parser.add_argument("--output", type=Path, default=ROOT / "artifacts" / "bs-greeks")
    args = parser.parse_args()
    if not (args.write_golden or args.cpp):
        parser.error("pass --write-golden or --cpp")
    if args.write_golden:
        write_fixture(ROOT / "testdata" / "bs_greeks_golden.csv", fixtures(GOLDEN))
    if args.cpp:
        return compare(args.cpp.resolve(), args.output)


if __name__ == "__main__":
    sys.exit(main())
