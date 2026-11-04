# Exact check of [B1] (1),(2) = [Bir00b] (3),(4), and of [Bir00b] (14).
# Gaussian rationals via Fraction pairs, so no rounding at all.
# A test on random systems, NOT a proof; the proofs are in Appendix A of the note.
from fractions import Fraction as F
import random

class G:
    def __init__(s, a, b=0): s.a, s.b = F(a), F(b)
    def __add__(s, o): o = o if isinstance(o, G) else G(o); return G(s.a+o.a, s.b+o.b)
    __radd__ = __add__
    def __sub__(s, o): o = o if isinstance(o, G) else G(o); return G(s.a-o.a, s.b-o.b)
    def __mul__(s, o): o = o if isinstance(o, G) else G(o); return G(s.a*o.a-s.b*o.b, s.a*o.b+s.b*o.a)
    __rmul__ = __mul__
    def __eq__(s, o): o = o if isinstance(o, G) else G(o); return s.a == o.a and s.b == o.b

def run(n):
    z = [G(1)] + [G(F(random.randint(-9, 9), random.randint(1, 5)), F(random.randint(-9, 9), random.randint(1, 5))) for _ in range(n-1)]
    S = [None] + [sum((zi_pow(zi, j) for zi in z), G(0)) for j in range(1, n+1)]
    # b_j: (z-z_2)...(z-z_n) = z^(n-1) + b_1 z^(n-2) + ... + b_(n-1)
    poly = [G(1)]
    for r in z[1:]:
        new = poly + [G(0)]
        for i in range(1, len(new)): new[i] = new[i] - r*poly[i-1]
        poly = new
    b = poly  # b[0] = 1
    ok = True
    for t in range(1, n):  # (3)
        lhs = sum((b[i]*S[t-i] for i in range(t)), G(0))
        rhs = sum((b[i] for i in range(t)), G(0)) - t*b[t]
        ok &= (lhs == rhs)
    lhs = sum((b[i]*S[n-i] for i in range(n)), G(0))  # (4)
    ok &= (lhs == sum((b[i] for i in range(n)), G(0)))
    c = [sum((b[i] for i in range(t+1)), G(0)) for t in range(n)]  # c_0 = 1
    for t in range(1, n):  # (14), valid for every t <= n-1
        lhs = sum((c[i]*S[t-i] for i in range(t)), G(0))
        rhs = 2*sum((c[i] for i in range(t)), G(0)) - t*c[t]
        ok &= (lhs == rhs)
    return ok

def zi_pow(x, j):
    r = G(1)
    for _ in range(j): r = r*x
    return r

random.seed(1)
res = [run(n) for n in range(2, 10) for _ in range(20)]
print("systems tested:", len(res), "all identities exact:", all(res))
