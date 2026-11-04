"""Black-Scholes prices, Greeks and implied volatility for European options, using only the standard library.

Rates and the dividend yield are continuously compounded. Vega and rho are per unit change (multiply by 0.01
for a one point move); theta is per year (divide by 365 for a calendar day).
"""
from math import erf, exp, log, pi, sqrt


def norm_cdf(x):
    return 0.5 * (1.0 + erf(x / sqrt(2.0)))


def norm_pdf(x):
    return exp(-0.5 * x * x) / sqrt(2.0 * pi)


def _d1_d2(S, K, T, sigma, r, q):
    if S <= 0 or K <= 0 or T <= 0 or sigma <= 0:
        raise ValueError("S, K, T and sigma must be positive")
    v = sigma * sqrt(T)
    d1 = (log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / v
    return d1, d1 - v


def price(S, K, T, sigma, r=0.0, q=0.0, call=True):
    d1, d2 = _d1_d2(S, K, T, sigma, r, q)
    if call:
        return S * exp(-q * T) * norm_cdf(d1) - K * exp(-r * T) * norm_cdf(d2)
    return K * exp(-r * T) * norm_cdf(-d2) - S * exp(-q * T) * norm_cdf(-d1)


def greeks(S, K, T, sigma, r=0.0, q=0.0, call=True):
    d1, d2 = _d1_d2(S, K, T, sigma, r, q)
    dq, dr = exp(-q * T), exp(-r * T)
    gamma = dq * norm_pdf(d1) / (S * sigma * sqrt(T))
    vega = S * dq * norm_pdf(d1) * sqrt(T)
    decay = -S * dq * norm_pdf(d1) * sigma / (2.0 * sqrt(T))
    if call:
        return {"delta": dq * norm_cdf(d1), "gamma": gamma, "vega": vega,
                "theta": decay - r * K * dr * norm_cdf(d2) + q * S * dq * norm_cdf(d1),
                "rho": K * T * dr * norm_cdf(d2)}
    return {"delta": -dq * norm_cdf(-d1), "gamma": gamma, "vega": vega,
            "theta": decay + r * K * dr * norm_cdf(-d2) - q * S * dq * norm_cdf(-d1),
            "rho": -K * T * dr * norm_cdf(-d2)}


def implied_vol(target, S, K, T, r=0.0, q=0.0, call=True, lo=1e-6, hi=5.0, tol=1e-12):
    """Volatility that reproduces a market price, by bisection (price is increasing in sigma)."""
    f = lambda s: price(S, K, T, s, r, q, call) - target
    if f(lo) > 0 or f(hi) < 0:
        raise ValueError("no volatility in [lo, hi] gives this price")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


if __name__ == "__main__":
    c = price(100, 100, 1, 0.2, 0.05)
    print(f"call {c:.4f}, put {price(100, 100, 1, 0.2, 0.05, call=False):.4f}, implied vol {implied_vol(c, 100, 100, 1, 0.05):.4f}")
