"""
check_phase2c_algebra.py  --  exact (sympy) checks of the algebra in Lemma 8
(explicit_constants.md, sections 2.12-2.16).

Working note, unrefereed. The method is A. Biro's (Indag. Math. 11 (2000) 343-358);
we only make his constants explicit.

What this checks: pure algebraic identities used in the written proofs.
What it does NOT check: the inequalities themselves (monotonicity, triangle
inequalities, Lemma 6, Lemma 7). Those are argued in the text.
Part F at the end is a float sanity test (random numbers), not a proof.
"""
import random
import sympy as sp

ok_all = True


def check(name, expr):
    global ok_all
    ok = sp.simplify(expr) == 0
    ok_all = ok_all and ok
    print(("OK   " if ok else "FAIL ") + name)


h, C, c, x, H, phi, mu, q, rho, r, t, beta0 = sp.symbols(
    "h C c x H phi mu q rho r t beta0", real=True)
F = lambda y: (2 * y**2 + 1) / (3 * y)

# A. (18'): one step of the h-recursion
check("A1 (18') (h C + c x)/(C + c) = h + (x - h) c/(C + c)",
      (h * C + c * x) / (C + c) - (h + (x - h) * c / (C + c)))
check("A2 induction step: (h C + x c)/(C+c) - H = ((h-H) C + (x-H) c)/(C+c)",
      (h * C + x * c) / (C + c) - H - ((h - H) * C + (x - H) * c) / (C + c))
check("A3 1/(sqrt2 (t+1)) - 1/(2 sqrt2 t) = (t-1)/(2 sqrt2 t (t+1))",
      1 / (sp.sqrt(2) * (t + 1)) - 1 / (2 * sp.sqrt(2) * t)
      - (t - 1) / (2 * sp.sqrt(2) * t * (t + 1)))

# B. (19'): Lemma 1(ii) with |w1| = 2h, |w2| = 3/2 + rho, cos = x
f7 = q**2 - sp.Rational(1, 4)
lhs = (2 * h)**2 + (sp.Rational(3, 2) + rho)**2 - 2 * (2 * h) * (sp.Rational(3, 2) + rho) * x
T = 3 * rho + rho**2 - 4 * h * rho * x
check("B1 q^2 - |w1-w2|^2 = 6 h x - (4 h^2 + 2 - f7 + T)",
      (q**2 - lhs) - (6 * h * x - (4 * h**2 + 2 - f7 + T)))
check("B2 case rho >= 0:  T + rho = rho^2 + 4 rho (1 - h x)",
      (T + rho) - (rho**2 + 4 * rho * (1 - h * x)))
Tneg = T.subs(rho, -r)
check("B3 case rho = -r < 0:  T + r = r^2 + r (4 h x - 2)",
      (Tneg + r) - (r**2 + r * (4 * h * x - 2)))
R = sp.symbols("R", real=True)
check("B4 (4h^2 + 2 - f7 - R)/(6h) = F(h) - (f7 + R)/(6h)",
      (4 * h**2 + 2 - f7 - R) / (6 * h) - (F(h) - (f7 + R) / (6 * h)))
check("B5 cos^2(pi/8) = (2 + sqrt2)/4", sp.cos(sp.pi / 8)**2 - (2 + sp.sqrt(2)) / 4)
check("B6 4 cos^2(pi/8) - 2 = sqrt2", 4 * sp.cos(sp.pi / 8)**2 - 2 - sp.sqrt(2))
check("B7 4 cos^2(asin(q)/2) - 2 = 2 sqrt(1 - q^2)",
      sp.simplify(4 * sp.cos(sp.asin(q) / 2)**2 - 2 - 2 * sp.sqrt(1 - q**2)))
check("B8 f7 = q^2 - 1/4 = (q - 1/2)(q + 1/2)", f7 - (q - sp.Rational(1, 2)) * (q + sp.Rational(1, 2)))

# C. (20) and the two cases
check("C1 F'(h) = (2h^2 - 1)/(3h^2)", sp.diff(F(h), h) - (2 * h**2 - 1) / (3 * h**2))
check("C2 (F(h) - h)' = -(1 + h^2)/(3h^2)", sp.diff(F(h) - h, h) + (1 + h**2) / (3 * h**2))
check("C3 F(H) - H = (1 - H^2)/(3H)", F(H) - H - (1 - H**2) / (3 * H))
check("C4 with H = 1 - phi: F(H) - H - 2 phi/3 = phi^2/(3 (1 - phi))",
      (F(1 - phi) - (1 - phi)) - 2 * phi / 3 - phi**2 / (3 * (1 - phi)))
check("C5 case A: (1 - phi/3 - mu phi) - (1 - phi) = phi (2/3 - mu)",
      (1 - phi / 3 - mu * phi) - (1 - phi) - phi * (sp.Rational(2, 3) - mu))
check("C6 case B: 2 phi/3 - mu phi = K7 phi with K7 = 2/3 - mu",
      2 * phi / 3 - mu * phi - (sp.Rational(2, 3) - mu) * phi)
check("C7 Lemma 8(ii): H' + 2(1 - H')/3 = 1 - (1 - H')/3",
      H + 2 * (1 - H) / 3 - (1 - (1 - H) / 3))

# D. constants
K7 = sp.Rational(2, 3) - mu
K9 = sp.sqrt(2) * (1 - sp.cos(beta0)) / K7
K10 = 1 / (1 + q)
K11 = sp.sqrt(2) * (1 + q) * (1 - sp.cos(beta0)) / K7
check("D1 K11 = K9 / K10", K11 - K9 / K10)
check("D2 K12 = 1/(1 + K2^(i)) = 1/(2 + q)", 1 / (1 + (1 + q)) - 1 / (2 + q))

# E. how the contradiction closes: phi >= K11/L and L' > L > 0 give phi L' > K11
Lsym, Lp = sp.symbols("L Lp", positive=True)
check("E1 (K11/L) L' - K11 = K11 (L' - L)/L", (K11 / Lsym) * Lp - K11 - K11 * (Lp - Lsym) / Lsym)

# F. float sanity test of (19') in abstract form (NOT a proof)
random.seed(1)
import math
worst = float("inf")
n_tested = 0
c0 = math.cos(math.pi / 8)
for _ in range(200000):
    hh = random.uniform(c0, 1.0)
    bb = random.uniform(0.0, math.pi / 8)
    rr = random.uniform(-0.5, 0.5)
    w1 = 2 * hh
    w2 = complex(1.5 + rr, 0) * complex(math.cos(bb), -math.sin(bb))
    d2 = abs(w1 - w2) ** 2
    rhs = F(hh) - ((d2 - 0.25) + abs(rr)) / (6 * hh)
    worst = min(worst, math.cos(bb) - rhs)
    n_tested += 1
okF = worst >= -1e-12
ok_all = ok_all and okF
print(("OK   " if okF else "FAIL ") +
      f"F  float sanity: cos(beta) >= F(h) - (|w1-w2|^2 - 1/4 + |rho|)/(6h) on {n_tested} random cases; "
      f"min slack = {worst:.3e}")

print("ALL OK" if ok_all else "SOME CHECK FAILED")
