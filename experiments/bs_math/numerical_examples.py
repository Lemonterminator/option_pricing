"""Reproduce the paper workbook's numbers and figures; no production API changes.

Python 3, NumPy, SciPy (existing kernels), Matplotlib and mpmath are required.
High precision comparisons use an independent complete complex CDF, not Phi_t.
"""
from pathlib import Path
import csv
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "artifacts" / "bs-math"
sys.path.insert(0, str(ROOT / "Python_prototype"))
import black_scholes as production
from bs_batch import evaluate_call

mp.mp.dps = 80
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": .25,
                     "figure.dpi": 140, "savefig.dpi": 200})


def Phi(x):
    return mp.erfc(-x / mp.sqrt(2)) / 2


def density(x):
    return mp.exp(-x*x/2) / mp.sqrt(2*mp.pi)


def price(args, kind="call", order=None):
    S, K, r, q, vol, T = args
    tau = vol * mp.sqrt(T)
    d1 = (mp.log(S/K) + (r-q)*T) / tau + tau/2
    d2 = d1-tau

    def cdf(z):
        if order is None:
            return Phi(z)
        a, b = mp.re(z), mp.im(z)
        # Independent derivative definition supports diagnostic P1 as well.
        return Phi(a) + sum(mp.diff(Phi, a, n)*(1j*b)**n/mp.factorial(n)
                            for n in range(1, order+1))
    A, B = S*mp.exp(-q*T), K*mp.exp(-r*T)
    if kind == "call":
        return A*cdf(d1)-B*cdf(d2)
    return B*cdf(-d2)-A*cdf(-d1)


def greeks(args, kind="call"):
    S, K, r, q, vol, T = args
    tau = vol * mp.sqrt(T)
    d1 = (mp.log(S/K)+(r-q)*T)/tau+tau/2
    d2 = d1-tau
    A, B = S*mp.exp(-q*T), K*mp.exp(-r*T)
    gamma = mp.exp(-q*T)*density(d1)/(S*tau)
    vega = A*density(d1)*mp.sqrt(T)
    if kind == "call":
        delta = mp.exp(-q*T)*Phi(d1)
        rho, psi = T*B*Phi(d2), -T*A*Phi(d1)
        theta = -A*density(d1)*vol/(2*mp.sqrt(T))-r*B*Phi(d2)+q*A*Phi(d1)
        dual_delta = -mp.exp(-r*T)*Phi(d2)
    else:
        delta = -mp.exp(-q*T)*Phi(-d1)
        rho, psi = -T*B*Phi(-d2), T*A*Phi(-d1)
        theta = -A*density(d1)*vol/(2*mp.sqrt(T))+r*B*Phi(-d2)-q*A*Phi(-d1)
        dual_delta = mp.exp(-r*T)*Phi(-d2)
    g = (r-q)/tau-d2/(2*T)
    p = d1*d2
    return dict(price=price(args, kind), delta=delta, gamma=gamma, vega=vega,
                theta=theta, rho=rho, psi=psi,
                vanna=-mp.exp(-q*T)*density(d1)*d2/vol,
                vomma=vega*p/vol, speed=-gamma/S*(1+d1/tau),
                zomma=gamma/vol*(p-1), charm=q*delta-mp.exp(-q*T)*density(d1)*g,
                color=gamma*(q+1/(2*T)+d1*g), veta=vega*(q-1/(2*T)+d1*g),
                ultima=vega/vol**2*(p*(p-1)-d1*d1-d2*d2),
                dual_delta=dual_delta, dual_gamma=mp.exp(-r*T)*density(d2)/(K*tau))


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def precise_args(values):
    return tuple(mp.mpf(str(x)) for x in values)


def references():
    scenarios = {
        "scenario1": precise_args((100, 100, 0, 0, .2, 1)),
        "scenario2": precise_args((100, 100, 0, 0, .01, 1))[:-1] + (mp.mpf(1)/365,),
        "carry_atm": precise_args((100, 100, .03, .01, .2, 1)),
        "carry_otm": precise_args((95, 100, .03, .01, .25, .5)),
    }
    refs, checks = [], []
    # Independent derivatives of the full high-precision price, teaching evidence.
    indices = {
        "delta": ((0, 1),), "gamma": ((0, 2),), "vega": ((4, 1),),
        "theta": ((5, 1),), "rho": ((2, 1),), "psi": ((3, 1),),
        "vanna": ((0, 1), (4, 1)), "vomma": ((4, 2),),
        "speed": ((0, 3),), "zomma": ((0, 2), (4, 1)),
        "charm": ((0, 1), (5, 1)), "color": ((0, 2), (5, 1)),
        "veta": ((4, 1), (5, 1)), "ultima": ((4, 3),),
        "dual_delta": ((1, 1),), "dual_gamma": ((1, 2),),
    }
    for name, args in scenarios.items():
        for kind in ("call", "put"):
            values = greeks(args, kind)
            refs.append(dict(scenario=name, kind=kind,
                             **{k: mp.nstr(v, 24) for k, v in values.items()}))
            for greek, specs in indices.items():
                orders = [0]*6
                for i, n in specs:
                    orders[i] = n
                numerical = mp.diff(lambda *xs: price(xs, kind), args, tuple(orders))
                if greek in ("theta", "charm", "color", "veta"):
                    numerical = -numerical
                checks.append(dict(scenario=name, kind=kind, greek=greek,
                                   formula=mp.nstr(values[greek], 24),
                                   full_price_derivative=mp.nstr(numerical, 24),
                                   absolute_difference=mp.nstr(abs(values[greek]-numerical), 8)))
            S, K, r, q, vol, T = args
            residual = values["theta"]+vol**2*S**2*values["gamma"]/2+(r-q)*S*values["delta"]-r*values["price"]
            checks.append(dict(scenario=name, kind=kind, greek="PDE_residual",
                               formula="0", full_price_derivative=mp.nstr(residual, 24),
                               absolute_difference=mp.nstr(abs(residual), 8)))
    write_csv(OUT / "greeks_reference.csv", refs)
    write_csv(OUT / "derivative_comparisons.csv", checks)
    # P1 CS-real inconsistency, using tiny but resolvable 80-digit steps.
    args = scenarios["scenario1"]
    truth = greeks(args)["gamma"]
    rows = []
    for exponent in (1, 2, 3, 4, 5, 6):
        h = mp.mpf(10)**(-exponent)
        shifted = (args[0]+1j*h,) + args[1:]
        for order in (1, 2, 3):
            estimate = -2*(mp.re(price(shifted, order=order))-price(args))/h**2
            rows.append(dict(h=mp.nstr(h, 8), cdf_order=order,
                             gamma_reference=mp.nstr(truth, 24),
                             cs_real_gamma=mp.nstr(estimate, 24),
                             estimate_over_gamma=mp.nstr(estimate/truth, 16)))
    write_csv(OUT / "p1_gamma_limit.csv", rows)
    maxdiff = max(mp.mpf(row["absolute_difference"]) for row in checks)
    print(f"80-digit formula/derivative comparisons: {len(checks)}, max absolute difference {mp.nstr(maxdiff, 5)}")
    for kind in ("call", "put"):
        print("carry_atm", kind, {k:mp.nstr(v, 13) for k,v in greeks(scenarios["carry_atm"], kind).items()})


def numeric_sweeps():
    scenarios = [("scenario1", (100.,100.,0.,0.,.2,1.)),
                 ("scenario2", (100.,100.,0.,0.,.01,1/365))]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    summary = []
    for rowidx, (name, args) in enumerate(scenarios):
        delta, gamma = float(production.delta_call(*args)), float(production.Gamma(*args))
        rows = []
        for hrel in np.logspace(-16, -4, 24):
            h = args[0]*hrel
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                fd = float(production.delta_forward(*args, h=h))
                cs = float(production.delta_cs(*args, h=h, order=3))
                gf = float(production.gamma_forward(*args, h=h))
                gr = float(production.gamma_cs_real(*args, h=h, order=3))
                g45 = float(production.gamma_45_deg(*args, h=h, order=3))
            rows.append(dict(h_rel=hrel, h=h, Delta_analytic=delta, Delta_fd=fd, Delta_cs=cs,
                             err_D_fd=abs(fd-delta), err_D_cs=abs(cs-delta), Gamma_analytic=gamma,
                             Gamma_fd=gf, Gamma_cs_real=gr, Gamma_cs_45=g45,
                             err_G_fd=abs(gf-gamma), err_G_cs_real=abs(gr-gamma), err_G_cs_45=abs(g45-gamma)))
        write_csv(OUT / f"{name}.csv", rows)
        for col, keys in enumerate((("err_D_fd", "err_D_cs"), ("err_G_fd", "err_G_cs_real", "err_G_cs_45"))):
            ax = axes[rowidx,col]
            for key in keys:
                vals = np.array([row[key] for row in rows])
                hs = np.array([row["h_rel"] for row in rows])
                mask = np.isfinite(vals) & (vals>0)
                ax.loglog(hs[mask], vals[mask], ".-", label=key)
                valid = vals[np.isfinite(vals)]
                if len(valid):
                    best = np.flatnonzero(np.isfinite(vals))[np.argmin(valid)]
                    summary.append(dict(scenario=name, method=key, cdf_order=3, total=len(vals),
                                        finite=len(valid), nonfinite=len(vals)-len(valid), zero_error=int(np.sum(valid==0)),
                                        max=float(np.max(valid)), median=float(np.median(valid)),
                                        p99=float(np.percentile(valid,99)), min=float(np.min(valid)),
                                        best_h_rel=float(hs[best])))
            ax.set_title(f"{name}: {'Delta' if col==0 else 'Gamma'}")
            ax.set_xlabel("relative spot step h/S")
            ax.set_ylabel("absolute error vs analytic float64")
            ax.legend(fontsize=8)
    fig.savefig(OUT / "assignment_errors.png")
    plt.close(fig)
    write_csv(OUT / "error_summary.csv", summary)


def figures():
    S = np.linspace(50,150,301)
    values = evaluate_call(S,100,0,0,.2,1)
    fig, axes = plt.subplots(1,3,figsize=(10,3), constrained_layout=True)
    for ax, y, title in zip(axes, (values.price, values.delta, values.gamma), ("Call price", "Delta = price slope", "Gamma = Delta slope")):
        ax.plot(S, y); ax.set_title(title); ax.set_xlabel("Spot S")
    fig.savefig(OUT / "price_delta_gamma.png"); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7,3.3), constrained_layout=True)
    S = np.linspace(99.5,100.5,1001)
    for vol,T in ((.01,1),(.01,1/30),(.01,1/365)):
        ax.plot(S, evaluate_call(S,100,0,0,vol,T).gamma, label=f"vol={vol}, T={T:.5g}")
    ax.set(xlabel="Spot S",ylabel="Gamma",title="Local width is approximately S * vol * sqrt(T)")
    ax.legend(); fig.savefig(OUT / "gamma_width.png"); plt.close(fig)
    T = np.geomspace(1/365,2,160)
    gamma = evaluate_call(100,100,0,0,.2,T).gamma
    vega = .2*T*100**2*gamma
    theta = -.5*.2**2*100**2*gamma
    fig, axes = plt.subplots(1,3,figsize=(10,3.3),constrained_layout=True)
    for ax,y,title in zip(axes,(gamma,vega,theta),("Gamma","Vega per 1.00 vol","Theta per year")):
        ax.semilogx(T,y); ax.set(xlabel="Remaining years T",title=title)
    fig.savefig(OUT / "greeks_term_structure.png"); plt.close(fig)
    args = precise_args((95,100,.03,.01,.25,.5))
    g = greeks(args)
    alpha = np.linspace(-1,1,101)
    exact, linear, quadratic = [], [], []
    for a in alpha:
        s,v = mp.mpf(str(a))*3, mp.mpf(str(a))*.02
        changed = (args[0]+s,) + args[1:4] + (args[4]+v,args[5])
        lin = g["delta"]*s+g["vega"]*v
        quad = lin+g["gamma"]*s*s/2+g["vanna"]*s*v+g["vomma"]*v*v/2
        exact.append(float(price(changed)-g["price"]))
        linear.append(float(lin)); quadratic.append(float(quad))
    fig, axes = plt.subplots(1,2,figsize=(9,3.4),constrained_layout=True)
    axes[0].plot(alpha,exact,label="exact reprice")
    axes[0].plot(alpha,linear,"--",label="Delta + Vega")
    axes[0].plot(alpha,quadratic,":",label="second order incl Vanna")
    axes[0].set(xlabel="shock scale: dS=3a, dvol=0.02a",ylabel="price change")
    axes[0].legend(fontsize=8)
    axes[1].plot(alpha,np.array(linear)-exact,label="linear error")
    axes[1].plot(alpha,np.array(quadratic)-exact,label="second-order error")
    axes[1].set(xlabel="shock scale a",ylabel="signed approximation error")
    axes[1].legend(fontsize=8)
    fig.savefig(OUT / "taylor_risk_approximation.png"); plt.close(fig)
    write_csv(OUT / "taylor_shocks.csv", [dict(a=a,exact=e,linear=l,quadratic=q) for a,e,l,q in zip(alpha,exact,linear,quadratic)])


def measure_example():
    """A hypothetical P/Q terminal-event comparison, not a market forecast."""
    OUT.mkdir(parents=True,exist_ok=True)
    args = precise_args((100,100,.03,.01,.20,1))
    S,K,r,q,vol,T = args
    d2q = (mp.log(S/K)+(r-q-vol**2/2)*T)/(vol*mp.sqrt(T))
    rows=[]
    for mu in map(mp.mpf,(".08",".03",".01")):
        d2p=d2q+(mu-r)*mp.sqrt(T)/vol
        g=greeks(args)
        row=dict(S=S,K=K,r=r,q=q,vol=vol,T=T,mu_total_return=mu,
                 d2_Q=d2q,d2_P=d2p,q_prob_terminal_itm_bs=Phi(d2q),
                 p_prob_terminal_itm_model=Phi(d2p),Phi_d1_Q=Phi(d2q+vol*mp.sqrt(T)),
                 delta=g["delta"],call_price=g["price"])
        rows.append({key:mp.nstr(value,30) for key,value in row.items()})
    write_csv(OUT/"measure_probabilities.csv",rows)
    mu=mp.mpf(".08")
    center_q=float((r-q-vol**2/2)*T)
    center_p=float((mu-q-vol**2/2)*T)
    width=float(vol*mp.sqrt(T)); threshold=float(mp.log(K/S))
    x=np.linspace(-.7,.8,901)
    fig,axes=plt.subplots(1,2,figsize=(9,3.3),constrained_layout=True,sharex=True,sharey=True)
    for ax,center,label,color,prob in zip(axes,(center_q,center_p),("Q: pricing weights","P: assumed total return 8%"),
                                        ("#276d83","#bd702d"),(Phi(d2q),Phi(d2q+(mu-r)*mp.sqrt(T)/vol))):
        y=np.exp(-.5*((x-center)/width)**2)/(width*np.sqrt(2*np.pi))
        ax.plot(x,y,color=color)
        ax.fill_between(x,y,where=x>threshold,color=color,alpha=.24)
        ax.axvline(threshold,color="black",linestyle="--",linewidth=1)
        ax.set(xlabel="Log price return log(S_expiry / S)",title=f"{label}\nPr(S_expiry > 100) = {float(prob):.4%}")
    axes[0].set_ylabel("Probability density")
    fig.savefig(OUT/"measure_distributions.png"); plt.close(fig)
    print("P/Q teaching probabilities: "+str(OUT/"measure_probabilities.csv"))


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    references()
    numeric_sweeps()
    figures()
    measure_example()
    print(f"Teaching artifacts: {OUT}")


if __name__ == "__main__":
    main()
