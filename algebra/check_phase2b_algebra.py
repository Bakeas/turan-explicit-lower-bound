"""Symbolic checks of the algebra used in explicit_constants.md, part 2
(Lemma 6, formula (14), Lemma 7 of Biro [Bir00b]). Exact, sympy.
Only the algebra between lines is checked here; the inequalities are argued in the note."""
import sympy as sp

ok_all = True
def check(name, expr):
    global ok_all
    e = sp.simplify(sp.expand(expr))
    good = (e == 0)
    ok_all &= good
    print(("OK  " if good else "FAIL"), name, "" if good else f"-> {e}")

s2 = sp.sqrt(2)
half = sp.Rational(1, 2)
q, B, Bt, t, k, e, e1, r, EA, Ec, f3 = sp.symbols('q B B_t t k e e1 r E_A E_c f3', real=True)
delta = q - half
lam = s2*delta  # Lemma 4 constant: lambda = sqrt2 (q - 1/2)

# (13): sum_{j=1}^t j|b_j| + sum_{j=1}^t B_{j-1} = t B_t, with B_j = 1 + |b_1| + ... + |b_j|
for T in range(1, 9):
    bb = sp.symbols(f'a1:{T+1}', nonnegative=True)  # a_j = |b_j|
    Bs = [1 + sum(bb[:j]) for j in range(T + 1)]     # Bs[j] = B_j
    check(f"(13) identity, t={T}", sum((j+1)*bb[j] for j in range(T)) + sum(Bs[j] for j in range(T)) - T*Bs[T])

# Lemma 6, t <= k-1:  A_t = (2/3) t B_t + E_A,  C_{t-1} = A_t/sqrt2 + E_c,  |c_t| = B_t/sqrt2 + e1
A_t = sp.Rational(2, 3)*t*Bt + EA
C = A_t/s2 + Ec
check("Lemma 6 (t<k): t|c_t| - (3/2)C = t e1 - (3/2)E_c - (3/(2 sqrt2)) E_A",
      t*(Bt/s2 + e1) - sp.Rational(3, 2)*C - (t*e1 - sp.Rational(3, 2)*Ec - 3/(2*s2)*EA))
check("f_4: (3/sqrt2) C - t B_t = (3/2) E_A + (3/sqrt2) E_c",
      3/s2*C - t*Bt - (sp.Rational(3, 2)*EA + 3/s2*Ec))
# bound assembly: (3/2) lam + (3/(2 sqrt2)) (2/3) f3 == (3 delta + f3)/sqrt2
check("Lemma 6 (t<k): (3/2)lam + (1/sqrt2) f3 = (3 delta + f3)/sqrt2",
      sp.Rational(3, 2)*lam + 3/(2*s2)*sp.Rational(2, 3)*f3 - (3*delta + f3)/s2)
# K_4 = 1/(sqrt2 (1 + K_2)), K_2 = 1 + q  =>  1/(sqrt2 K_4) = 2 + q
K2 = 1 + q
K4 = 1/(s2*(1 + K2))
check("1/(sqrt2 K_4) = 2 + q", 1/(s2*K4) - (2 + q))

# Lemma 6, t = k:  A_{k-1} = (2/3)(k-1)B + E_A,  A_k = A_{k-1} + B,  C_{k-1} = A_k/sqrt2 + E_c,
#                  k|c_k| = k B/sqrt2 + k e + r      (B = B_{k-1})
A_km1 = sp.Rational(2, 3)*(k - 1)*B + EA
A_k = A_km1 + B
Ck = A_k/s2 + Ec
check("Lemma 6 (t=k): k|c_k| - (3/2)C_{k-1} = (k e - (3/2)E_c) + r - B/(2 sqrt2) - (3/(2 sqrt2)) E_A",
      (k*B/s2 + k*e + r) - sp.Rational(3, 2)*Ck - ((k*e - sp.Rational(3, 2)*Ec) + r - B/(2*s2) - 3/(2*s2)*EA))
check("Lemma 6 (t=k): A_k = (2/3) k B + B/3 + E_A", A_k - (sp.Rational(2, 3)*k*B + B/3 + EA))
# bound assembly at t=k: (3/2)k lam + K2 + 1/(2 sqrt2) + (1/sqrt2)(k-1) f3  <=  k[(3 delta + f3)/sqrt2 + (K2 + 1/(2sqrt2))/k]
lhs = sp.Rational(3, 2)*k*lam + K2 + 1/(2*s2) + (k - 1)*f3/s2
rhs = k*((3*delta + f3)/s2 + (K2 + 1/(2*s2))/k)
check("Lemma 6 (t=k): rhs - lhs = f3/sqrt2 (>= 0)", rhs - lhs - f3/s2)
check("sqrt2 (K_2 + 1/(2 sqrt2)) = sqrt2 (1+q) + 1/2", s2*(K2 + 1/(2*s2)) - (s2*(1 + q) + half))

# (14): L_t = sum_{j=0}^{t-1} c_j S_{t-j},  R_t = 2 D_{t-1} - t c_t,  c_j = 1 + b_1 + ... + b_j
for T in range(2, 9):
    b = [sp.Integer(1)] + list(sp.symbols(f'b1:{T+1}'))
    S = [None] + list(sp.symbols(f'S1:{T+1}'))
    c = [sum(b[:j+1]) for j in range(T + 1)]
    L = lambda tt: sum(c[j]*S[tt - j] for j in range(tt))
    R = lambda tt: 2*sum(c[:tt]) - tt*c[tt]
    check(f"(14) step L_t - L_(t-1) = sum b_j S_(t-j), t={T}", L(T) - L(T-1) - sum(b[j]*S[T - j] for j in range(T)))
    check(f"(14) step R_t - R_(t-1) = c_(t-1) - t b_t,   t={T}", R(T) - R(T-1) - (c[T-1] - T*b[T]))
b1, S1 = sp.symbols('b1 S1')
check("(14) base t=1: R_1 = 2 - c_1 = 1 - b_1", 2 - (1 + b1) - (1 - b1))

# Lemma 7: 2 sin(beta) cos(beta) = sin(2 beta), and sin(2*pi/8) = 1/sqrt2
bt = sp.symbols('beta', real=True)
check("2 sin b cos b = sin 2b", 2*sp.sin(bt)*sp.cos(bt) - sp.sin(2*bt))
check("sin(pi/4) = 1/sqrt2", sp.sin(sp.pi/4) - 1/s2)

print("\nALL OK" if ok_all else "\nSOME CHECK FAILED")
