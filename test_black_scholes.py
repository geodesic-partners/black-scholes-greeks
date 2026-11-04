"""Run with: python -m pytest  (or python test_black_scholes.py)."""
import itertools
from math import exp

from black_scholes import greeks, implied_vol, price

CASES = list(itertools.product((80.0, 100.0, 125.0), (90.0, 110.0), (0.25, 2.0), (0.1, 0.45), (0.0, 0.05), (0.0, 0.03)))


def test_textbook_value():
    # S = K = 100, one year, 20% vol, 5% rate: the call is worth 10.4506 to four places.
    assert abs(price(100, 100, 1, 0.2, 0.05) - 10.450583572185565) < 1e-9


def test_put_call_parity():
    for S, K, T, s, r, q in CASES:
        lhs = price(S, K, T, s, r, q) - price(S, K, T, s, r, q, call=False)
        assert abs(lhs - (S * exp(-q * T) - K * exp(-r * T))) < 1e-9


def test_greeks_match_finite_differences():
    h = 1e-4
    for (S, K, T, s, r, q), call in itertools.product(CASES, (True, False)):
        p = lambda S=S, T=T, s=s, r=r: price(S, K, T, s, r, q, call)
        g = greeks(S, K, T, s, r, q, call)
        fd = {"delta": (p(S=S + h) - p(S=S - h)) / (2 * h),
              "gamma": (p(S=S + h) - 2 * p() + p(S=S - h)) / h ** 2,
              "vega": (p(s=s + h) - p(s=s - h)) / (2 * h),
              "theta": -(p(T=T + h) - p(T=T - h)) / (2 * h),
              "rho": (p(r=r + h) - p(r=r - h)) / (2 * h)}
        for k in fd:
            assert abs(g[k] - fd[k]) < 1e-4 * max(1.0, abs(fd[k])), (k, S, K, T, s, r, q, call)


def test_implied_vol_round_trip():
    # Deep in or out of the money with little time left, the price hardly moves with volatility (vega near zero),
    # so the volatility cannot be read back from it to any useful precision. Those cases are skipped.
    for (S, K, T, s, r, q), call in itertools.product(CASES, (True, False)):
        if greeks(S, K, T, s, r, q, call)["vega"] < 0.01:
            continue
        v = implied_vol(price(S, K, T, s, r, q, call), S, K, T, r, q, call)
        assert abs(v - s) < 1e-6, (S, K, T, s, r, q, call)


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("passed", name)
