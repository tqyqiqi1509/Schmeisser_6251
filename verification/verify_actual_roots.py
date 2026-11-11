#!/usr/bin/env python3
"""Exact rational checks for the actual-zero estimates in Section 3.

No interval or numerical library is used. Acceptance uses Fraction only.
Logarithms are enclosed by the atanh series after exact range reduction;
negative exponentials are bounded by reciprocals of positive Taylor sums.
"""
from fractions import Fraction as F
from math import factorial
import json


def log_bounds(x, terms=48):
    x = F(x)
    assert x >= 1
    k = 0
    while x >= 2:
        x /= 2
        k += 1

    def unit(y):
        t = (y-1)/(y+1)
        lo = 2*sum((t**(2*j+1)/F(2*j+1) for j in range(terms)), F(0))
        err = 2*t**(2*terms+1)/(F(2*terms+1)*(1-t*t))
        return lo, lo+err

    lo2, hi2 = unit(F(2))
    lo, hi = unit(x)
    return lo+k*lo2, hi+k*hi2


def exp_partial(x, degree):
    x = F(x)
    assert x >= 0
    term = total = F(1)
    for j in range(1, degree+1):
        term *= x/F(j)
        total += term
    return total


def exp_negative_upper(x):
    assert x > 0
    return 1/exp_partial(x, 160)


def run():
    m = F(6250)
    rows = [
        (F(0), F(1,4), F(659,1000)),
        (F(1,4), F(3,8), F(396,1000)),
        (F(3,8), F(1,2), F(346,1000)),
        (F(1,2), F(5,8), F(192,1000)),
        (F(5,8), F(3,4), F(54,1000)),
        (F(3,4), F(7,8), F(11,1000)),
        (F(7,8), F(19,20), F(6,1000)),
    ]
    diagnostics = []
    for l,h,bound in rows:
        dm, dp = max(F(1,10), 1-h*h), 1-l*l
        lm, lp = (m+2)*dm/2, (m+2)*dp/2
        logm, _ = log_bounds(lm)
        logp, _ = log_bounds(lp)
        dec = 7*lp/(m*logp)
        if l:
            dec = min(dec, 28/(m*l*l))
        inc = F(5,3)*h*h*m*exp_negative_upper((m-1)/(m+2)*logm/dp)
        end = exp_negative_upper(m/(m+2)*logm/dp)
        total = dec+inc+end
        assert total < bound < 1-h
        diagnostics.append({"interval":[str(l),str(h)],"rational_upper":str(bound),
                            "display_only_upper":float(total)})
    assert rows[0][0] == 0
    assert all(rows[j][1] == rows[j+1][0] for j in range(len(rows)-1))
    assert rows[-1][1]**2 > F(9,10)

    # Interior-to-boundary transition used in the energy proof.
    t = F(700)
    lo35, _ = log_bounds(t/20)
    _, hi700 = log_bounds(t)
    assert (t-3)/8 > 4*hi700
    assert 10*(1-3/t)*lo35 > 4*hi700
    assert 9*t/160 > 4*hi700
    assert F(9,10) > F(94,100)**2
    assert F(10240,243)*t*t*(t+1)/((t-1)**3*F(94,100)**3) < 60
    assert F(10,3)*(t+1)/t < 7
    assert 61/t < F(94,200)
    assert 1/(2*t**3) < F(1,4)
    assert F(16)**2/F(2)**17 < F(1,16)

    # Constants in the variance proof. Subsequent degree monotonicity is
    # proved by differentiation in the paper.
    n = F(2048)
    m = n-1
    assert log_bounds(F(2))[0] > F(56,81)
    assert (1-1/n)*(22*F(56,81)-1) > 14
    assert m*m/(m-1)**2 < F(101,100)
    assert (m+1)*m*m/(m-1)**3 < F(101,100)
    assert F(6619136,25515) < 300
    assert F(51712,135) < 400
    assert exp_partial(F(7161,2047),7) > 32
    assert exp_partial(F(63,16),12) > 48
    assert F(25,192)+F(27,56)+F(35,128)+F(1,48) == F(2437,2688) < 1
    assert F(2)**17 > 32*n
    assert F(5,2)**19 > 6*n*n
    assert F(2)**47 > 384*n**3
    assert exp_partial(F(1),4) > F(5,2)
    e_upper = exp_partial(F(1),12) + F(1,factorial(13))/(1-F(1,14))
    assert e_upper < F(25,9)
    return {"status":"PASS","arithmetic":"exact rational",
            "energy_degree_floor":6251,"variance_degree_floor":2048,
            "energy_intervals":diagnostics}


if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
