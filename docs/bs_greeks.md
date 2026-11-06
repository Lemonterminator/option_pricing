# European Vega, Theta and Rho

The scalar and batch APIs take `S, K, r, q, sigma, T`. Spot and strike must be positive, volatility and time to expiry must be nonnegative, and all inputs must be finite. Rates are continuously compounded, time is in years, and outputs are per underlying unit. The APIs cover European Call and Put sensitivities with constant inputs and continuous yield.

| Field | Definition | Unit |
| --- | --- | --- |
| `vega` | Partial derivative with respect to `sigma` | Per unit volatility. Multiply by 0.01 for one volatility point. |
| `theta` | Calendar time derivative, the negative of the derivative with respect to `T` | Per year. Divide by 365 for a local calendar day approximation under a 365-day convention. |
| `rho` | Partial derivative with respect to `r`, holding spot and dividend yield fixed | Per unit rate. Multiply by 0.01 for one percentage point or 0.0001 for one basis point. |

## Calling the code

Run the teaching prototype from the repository root with `python Python_prototype/greeks_prototype.py`. To call it or the vectorized API, add `Python_prototype/` to the import path as for the existing price API.

```python
from greeks_prototype import greeks_scalar
from bs_greeks import evaluate_greeks

call = greeks_scalar(100, 100, .03, .01, .20, 1)
puts = evaluate_greeks([90, 100, 110], 100, .03, .01, .20, 1, kind="put")
print(call.theta, puts.rho)
```

`GreeksResult` contains three NumPy arrays with the broadcast shape, including `()` for scalar inputs and empty arrays when the broadcast is empty. The scalar prototype is a self-contained NumPy/SciPy program.

The C++17 implementation is in `black_scholes_core` and declared in `BS_cpp/bs_batch.h`. The batch overload takes the existing `CallInputs` shape, with six equal-length vectors. The result fields are `std::vector<double>`.

```cpp
#include "bs_batch.h"
using namespace black_scholes;

auto call = evaluate_greeks(100, 100, .03, .01, .20, 1);
// inputs holds six equal-length vectors in the usual parameter order.
auto puts = evaluate_greeks(inputs, OptionKind::Put, BatchBackend::Serial);
```

`BatchBackend::ParUnseq` is also available. Threading depends on the standard library and build. Price, Delta and Gamma remain available via `evaluate_call`. The Put API here returns sensitivities only.

## Boundaries

Let $A=S e^{-qT}$, $B=K e^{-rT}$, and let $\eta$ be +1 for Call and -1 for Put. Let $g=\eta(qA-rB)$.

| Case | Vega | Theta | Rho |
| --- | --- | --- | --- |
| $\sigma=0$, $T>0$, deterministic ITM | 0 | $g$ | $\eta TB$ |
| $\sigma=0$, $T>0$, deterministic OTM | 0 | 0 | 0 |
| $\sigma=0$, $T>0$, forward ATM | $A\sqrt{T/(2\pi)}$, the right derivative | $\max(g,0)$, the calendar right derivative | $\eta TB/2$, midpoint convention at the rate kink |
| $T=0$, ITM | 0 | $g$ evaluated at $T=0$ | 0 |
| $T=0$, OTM | 0 | 0 | 0 |
| $T=0$, $S=K$, $\sigma=0$ | 0 | $\min(g,0)$ evaluated at $T=0$ | 0 |
| $T=0$, $S=K$, $\sigma>0$ | 0 | $-\infty$, the expiry limit | 0 |

Forward ATM at zero volatility is detected by a computed zero in $m=\log(S/K)+(r-q)T$. Values near a kink can fall on either side if the sign of $m$ is below floating point resolution. At expiry, Theta is defined by the limit from positive remaining maturity. At the zero-volatility rate kink the two-sided Rho does not exist. Its midpoint is also the positive-volatility limit when the other inputs are fixed at forward ATM.

Python raises `ValueError` for invalid domain inputs or kind and `OverflowError` when a calculation exceeds float64 range. C++ uses `std::invalid_argument` and `std::overflow_error`. The explicitly defined negative infinity at ATM expiry is returned as a value. A finite mathematical value can still underflow when rounded to double.

## Floating point evaluation

The kernel calculates $d_1$ and $d_2$ from log forward moneyness and combines discount, scale and normal tail in log space. The C++ path reuses the existing `Phi_real` on its usual range and uses the log of a Mills expansion below -35. This lets a sufficiently large notional retain a tail that would underflow if the CDF were rounded first. Python uses SciPy `log_ndtr`. Probit inverts the CDF and is not needed by the pricing derivatives.

For a narrow normal interval, calculate $\Phi(d_1)-\Phi(d_2)$ from the integral of the density instead of subtracting the CDF values. The symmetric expansion uses the midpoint $c=(d_1+d_2)/2$, interval width $\tau=\sigma\sqrt T$, and the correction $((c\tau)^2-\tau^2)/24$. The code takes this branch for $\tau(1+|c|)\le 10^{-3}$. This preserves small carry terms in Theta when $r$ and $q$ are nearly equal. See the [paper lesson](bs_math/06_first_order_greeks.md) for the derivation.

## Reproducing the checks

```powershell
python -m pip install -r Python_prototype/requirements.txt
python -m pip install mpmath QuantLib==1.44
python -m unittest discover -s tests -v
cmake -S BS_cpp -B out/greeks-cpp -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build out/greeks-cpp
ctest --test-dir out/greeks-cpp --output-on-failure
python -X utf8 experiments/bs_math/greeks_accuracy.py --cpp out/greeks-cpp/bs_greeks_contract_test.exe
```

Use `--write-golden` on the accuracy script to regenerate the [shared fixture](../testdata/bs_greeks_golden.csv). `--output` selects a different artifact directory. See [the measured accuracy report](bs_greeks_accuracy.md) for the cases, error scales and results.
