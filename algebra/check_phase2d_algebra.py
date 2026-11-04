"""
check_phase2d_algebra.py -- exact (sympy) checks of the algebraic identities used in
explicit_constants.md, part 4 (Section 4 of Biro's paper, proof of the Theorem; our 2.20-2.27).

Only identities are checked here. The inequalities (triangle inequality, Cauchy-Schwarz,
monotonicity) are argued in the text. Part F is a float sanity test (NOT a proof).

Method and ideas: A. Biro, Indag. Math. (N.S.) 11 (2000), 343-358.
Prepared by V. Bakeas with Claude (AI). Working note, unrefereed.
"""
import random
import sympy as sp

ok = True


def check(name, expr):
    global ok
    e = sp.simplify(sp.expand(expr))
    if e == 0:
        print(f"OK    {name}")
    else:
        ok = False
        print(f"FAIL  {name}: residual {e}")


s2 = sp.sqrt(2)
h, rho, psi, C, q, g, gam, f7, x, y = sp.symbols('h rho psi C q g gamma f7 x y', real=True)
I = sp.I

# ---------- A. (38): f_11 ----------
r = sp.Rational(3, 2) + rho
E = sp.exp(I * psi)
# A1: (3/2)h - r e^{i psi} = (3/2)(h-1) + (3/2)(1-e^{i psi}) - rho e^{i psi}
check("A1 (38) splitting identity",
      (sp.Rational(3, 2) * h - r * E)
      - (sp.Rational(3, 2) * (h - 1) + sp.Rational(3, 2) * (1 - E) - rho * E))
# A2: |1 - e^{i psi}|^2 = 2(1 - cos psi)
check("A2 |1-e^{i psi}|^2 = 2(1-cos psi)",
      sp.expand_complex((1 - E) * sp.conjugate(1 - E)) - 2 * (1 - sp.cos(psi)))
# A3: (14) at t=k minus D/2: (2D - k c_k) - D/2 = (3/2)D - k c_k
D, kc = sp.symbols('D kc')
check("A3 2D - kc_k - D/2 = (3/2)D - kc_k", (2 * D - kc - D / 2) - (sp.Rational(3, 2) * D - kc))
# A4: factorisation with D = h C u, k c_k = r C u e^{i psi}
u = sp.symbols('u')
check("A4 (3/2)D - k c_k = C u ((3/2)h - r e^{i psi})",
      (sp.Rational(3, 2) * h * C * u - r * C * u * E) - C * u * (sp.Rational(3, 2) * h - r * E))

# ---------- B. (40)-(41): f_12, f_13 ----------
check("B1 1/h - 1 = (1-h)/h", (1 / h - 1) - (1 - h) / h)
f12 = g / (1 - g)
check("B2 f12 at h = 1-g", (1 / (1 - g) - 1) - f12)
check("B3 f12(1+f12) = g/(1-g)^2", f12 * (1 + f12) - g / (1 - g) ** 2)
# B4: Lemma 1(v): ||w| - w|^2 = 2(|w| - Re w)|w|, w = x + i y, R := |w|
R = sp.symbols('R', positive=True)
w = x + I * y
lhs = sp.expand((R - w) * (R - sp.conjugate(w)))
rhs = 2 * (R - x) * R
check("B4 Lemma 1(v) using x^2 + y^2 = R^2", sp.expand(lhs - rhs).subs(x**2, R**2 - y**2))
# B5: dependence on h: sqrt(2(1/h-1)/h) = sqrt(2(1-h))/h  (squares)
check("B5 2(1/h-1)(1/h) = 2(1-h)/h^2", 2 * (1 / h - 1) * (1 / h) - 2 * (1 - h) / h**2)

# ---------- C. (42), H_c: f_14, f_15 ----------
TA, TH, cc, f12s, f14s = sp.symbols('T_all T_H c f12 f14', real=True)
# C1: q*(T_all - T_H) + c*T_H = q*T_all - (q - c)*T_H
check("C1 H_c rearrangement", q * (TA - TH) + cc * TH - (q * TA - (q - cc) * TH))
# C2: q(1 + f12) - 1/2 + f14 = (q - 1/2) + q f12 + f14
check("C2 f15 assembly", q * (1 + f12s) - sp.Rational(1, 2) + f14s
      - ((q - sp.Rational(1, 2)) + q * f12s + f14s))
# C3: with c = 1/2 - gamma: 1/2 - c = gamma; q - c = (q - 1/2) + gamma
check("C3a 1/2 - c = gamma", (sp.Rational(1, 2) - (sp.Rational(1, 2) - gam)) - gam)
check("C3b q - c = (q-1/2) + gamma", (q - (sp.Rational(1, 2) - gam)) - ((q - sp.Rational(1, 2)) + gam))

# ---------- D. f_18, f_19, (44) ----------
Sx, Sy = sp.symbols('S_x S_y', real=True)
S = Sx + I * Sy
# D1: |S - 1/2|^2 = |S|^2 - Re S + 1/4
check("D1 |S-1/2|^2 = |S|^2 - Re S + 1/4",
      sp.expand_complex((S - sp.Rational(1, 2)) * sp.conjugate(S - sp.Rational(1, 2)))
      - ((Sx**2 + Sy**2) - Sx + sp.Rational(1, 4)))
# D2: q^2 - (1/2 - gamma) + 1/4 = (q^2 - 1/4) + gamma
check("D2 q^2 - c + 1/4 = f7 + gamma",
      q**2 - (sp.Rational(1, 2) - gam) + sp.Rational(1, 4) - ((q**2 - sp.Rational(1, 4)) + gam))
# D3: for k = 2..7: sum_{j<k} b_j S_{k-j} - c_{k-1}/2 = sum_j b_j (S_{k-j} - 1/2), b_0 = 1
for k in range(2, 8):
    b = [sp.Integer(1)] + list(sp.symbols(f'b1:{k}'))
    Ss = sp.symbols(f'S1:{k + 1}')
    ck1 = sum(b)
    lhs = sum(b[j] * Ss[k - j - 1] for j in range(k)) - ck1 / 2
    rhs = sum(b[j] * (Ss[k - j - 1] - sp.Rational(1, 2)) for j in range(k))
    check(f"D3 splitting identity, k={k}", lhs - rhs)
# D4: (44): (c - kb) - c/2 = c/2 - kb
cS, kb = sp.symbols('c kb')
check("D4 (44) rearrangement", (cS - kb) - cS / 2 - (cS / 2 - kb))

# ---------- E. final step ----------
check("E1 1/(1/sqrt2 - 1/2) = 2 + 2 sqrt2", 1 / (1 / s2 - sp.Rational(1, 2)) - (2 + 2 * s2))
check("E2 (1/sqrt2)(1/sqrt2 - 1/2) = (2 - sqrt2)/4",
      (1 / s2) * (1 / s2 - sp.Rational(1, 2)) - (2 - s2) / 4)
a_, b_ = sp.symbols('a b', real=True)
check("E3 2(a^2+b^2) - (a+b)^2 = (a-b)^2", 2 * (a_**2 + b_**2) - (a_ + b_)**2 - (a_ - b_)**2)
check("E4 sqrt2 * 1/(2 sqrt2) = 1/2", s2 * (1 / (2 * s2)) - sp.Rational(1, 2))
check("E5 sqrt2*f19 < 1/(2 sqrt2) <=> f19 < 1/4: 1/(2 sqrt2)/sqrt2 = 1/4",
      (1 / (2 * s2)) / s2 - sp.Rational(1, 4))
# E6: k b_k / c = 1/2 + e with e = (k b_k - c/2)/c
check("E6 kb/c = 1/2 + (kb - c/2)/c", kb / cS - (sp.Rational(1, 2) + (kb - cS / 2) / cS))
# E7: (2 - sqrt2)/4 < 1/4
check("E7 1/4 - (2 - sqrt2)/4 = (sqrt2 - 1)/4 > 0", (sp.Rational(1, 4) - (2 - s2) / 4) - (s2 - 1) / 4)

# ---------- G. cancellation-free forms for the Phase 3 checker ----------
# x := pi/4 - alpha, d := q - 1/2
xx, d = sp.symbols('x d', real=True)
al = sp.pi / 4 - xx
qq = sp.Rational(1, 2) + d
def trig0(e):
    return sp.simplify(sp.expand_trig(sp.expand(e)).rewrite(sp.exp).expand()).rewrite(sp.cos)
check("G1 sin a = (cos x - sin x)/sqrt2", trig0(sp.sin(al) - (sp.cos(xx) - sp.sin(xx)) / s2))
check("G2 cos a - 1/sqrt2 = (sin x - 2 sin^2(x/2))/sqrt2",
      trig0(sp.cos(al) - 1 / s2 - (sp.sin(xx) - 2 * sp.sin(xx / 2)**2) / s2))
check("G3 q^2 - sin^2 a/2 = d + d^2 + sin(2x)/4",
      trig0(qq**2 - sp.sin(al)**2 / 2 - (d + d**2 + sp.sin(2 * xx) / 4)))
check("G4 q - sin a/sqrt2 = d + sin^2(x/2) + sin(x)/2",
      trig0(qq - sp.sin(al) / s2 - (d + sp.sin(xx / 2)**2 + sp.sin(xx) / 2)))
check("G5 f7 = q^2 - 1/4 = d(1 + d)", qq**2 - sp.Rational(1, 4) - d * (1 + d))
check("G6 M(q) = (sqrt2-1)/(2d): (1 - 1/sqrt2)/(sqrt2 d) - (sqrt2-1)/(2d)",
      (1 - 1 / s2) / (s2 * d) - (s2 - 1) / (2 * d))

# ---------- F. float sanity tests (NOT proofs) ----------
import cmath, math
random.seed(20261103)
worst41 = float('inf')
for _ in range(100000):
    m = random.randint(1, 12)
    t = [complex(random.gauss(0, 1), random.gauss(0, 1)) for _ in range(m)]
    s = sum(t)
    if abs(s) < 1e-3:
        continue
    t = [tj / s for tj in t]                      # now sum t_j = 1
    A = sum(abs(tj) for tj in t)
    lhs = sum(abs(abs(tj) - tj) for tj in t)
    rhs = math.sqrt(2 * (A - 1) * A)
    worst41 = min(worst41, rhs - lhs)
print(f"F1 float test of (41) bound on 100000 random cases: min slack {worst41:.3e} (should be >= -1e-12)")
if worst41 < -1e-12:
    ok = False
worst_cone = float('inf')
for _ in range(200000):
    r0 = random.random() * (1 / (2 * math.sqrt(2))) * 0.999999
    e = r0 * cmath.exp(1j * random.uniform(-math.pi, math.pi))
    wv = 0.5 + e
    worst_cone = min(worst_cone, math.pi / 4 - abs(cmath.phase(wv)))
print(f"F2 float test: |e| < 1/(2 sqrt2) => |arg(1/2+e)| < pi/4: min slack {worst_cone:.3e} (should be > 0)")
if worst_cone <= 0:
    ok = False
worst38 = float('inf')
for _ in range(200000):
    hh = random.uniform(0.5, 1)
    ps = random.uniform(-math.pi, math.pi)
    rr = random.uniform(-1, 1)
    lhs = abs(1.5 * hh - (1.5 + rr) * cmath.exp(1j * ps))
    rhs = 1.5 * (1 - hh) + 1.5 * math.sqrt(2 * (1 - math.cos(ps))) + abs(rr)
    worst38 = min(worst38, rhs - lhs)
print(f"F3 float test of (38) bound on 200000 random cases: min slack {worst38:.3e} (should be >= -1e-12)")
if worst38 < -1e-12:
    ok = False

print("ALL OK" if ok else "SOME CHECKS FAILED")
raise SystemExit(0 if ok else 1)
