# Black-Scholes prices, Greeks and implied volatility in Python

This is one file, `black_scholes.py`, which prices European calls and puts under Black-Scholes and gives you all five Greeks. You can also hand it a market price and it'll back out the implied volatility. It needs nothing beyond Python's standard library, so you can drop it into a notebook or an interview take-home without installing anything.

```python
from black_scholes import price, greeks, implied_vol

c = price(100, 100, 1, 0.2, r=0.05)          # 10.4506
greeks(100, 100, 1, 0.2, r=0.05)["delta"]     # 0.6368
implied_vol(c, 100, 100, 1, r=0.05)           # 0.2
```

Rates and the dividend yield `q` are continuously compounded. I've left vega and rho per unit change and theta per year, since that's what the formulas give you. Multiply vega and rho by 0.01 if you want a one point move, and divide theta by 365 if you'd like it per calendar day.

## Tests

`python -m pytest` (or just `python test_black_scholes.py`) checks four things over 96 combinations of the six inputs:

- the textbook value of 10.4506 for an at-the-money one-year call at 20% volatility and a 5% rate
- put-call parity to within 1e-9
- every Greek against a finite difference of the price, for calls and puts
- that implied volatility recovers the volatility a price was made with

The last test skips the 16 cases where vega is below 0.01. When an option is deep in or out of the money with only three months left, its price barely moves with volatility, and you can't read the volatility back from it with any real precision whatever method you use. I'd rather the test said so openly than pass by luck.

## Why implied volatility uses bisection

Newton's method on vega converges faster, but it can overshoot badly when vega is small. Bisection takes about 43 steps here and never leaves the bracket, which is a fair trade for a reference implementation.

If you're preparing for quant interviews, there are stochastic calculus and Black-Scholes questions with worked answers at https://geodesicpartners.co.uk/quant-questions-stochastic-calculus
