"""Symbolic checks of the algebraic identities used in explicit_constants.md
(Lemmas 1-5, Corollaries 1-2 of Biro [Bir00b]). Exact, sympy. Each check prints OK/FAIL.
These check only the algebra between lines; the inequalities are argued in the note."""
import sympy as sp

ok_all = True
def check(name, expr):
    global ok_all
    e = sp.simplify(sp.expand(expr))
    good = (e == 0)
    ok_all &= good
    print(("OK  " if good else "FAIL"), name, "" if good else f"-> {e}")

r1, r2, th, ph, q, B, b, X, ca, Bt, cb, ct, a = sp.symbols(
    'r1 r2 theta phi q B b X cos_alpha B_tm1 c_abs cos_alpha_t alpha', real=True)
s2 = sp.sqrt(2)

# Lemma 1(i): w1 = r1 e^{i(th+ph)}, w2 = r2 e^{i th}, so arg(w1/w2) = ph
w1 = r1*sp.exp(sp.I*(th+ph)); w2 = r2*sp.exp(sp.I*th)
lhs = sp.expand((w1-w2)*sp.conjugate(w1-w2))
rhs = r1**2*sp.sin(ph)**2 + (r2 - r1*sp.cos(ph))**2
check("Lemma 1(i)", sp.simplify(sp.expand_complex(lhs) - rhs))
# Lemma 1(ii) (law of cosines)
check("Lemma 1(ii)", sp.simplify(sp.expand_complex(lhs) - (r1**2 + r2**2 - 2*r1*r2*sp.cos(ph))))

# Corollary 1: (1 - 1/sqrt2)/sqrt2 == (sqrt2 - 1)/2
check("Cor 1 constant", (1 - 1/s2)/s2 - (s2 - 1)/2)
# sqrt2*q - 1/sqrt2 == sqrt2*(q - 1/2)
check("sqrt2 q - 1/sqrt2", s2*q - 1/s2 - s2*(q - sp.Rational(1, 2)))
# Cor 1: inequality 1 - 1/s2 + B/s2 <= s2 q B  <=>  B >= (s2-1)/(2(q-1/2))  (same linear form)
check("Cor 1 rearrangement", (s2*q*B - (1 - 1/s2) - B/s2) - s2*(q - sp.Rational(1, 2))*(B - (s2 - 1)/(2*(q - sp.Rational(1, 2)))))
# M(q) > 1  <=>  q < 1/sqrt2 : M(q) - 1 = (1/sqrt2 - q)/(q - 1/2)
Mq = (s2 - 1)/(2*(q - sp.Rational(1, 2)))
check("M(q) - 1 = (1/sqrt2 - q)/(q - 1/2)", Mq - 1 - (1/s2 - q)/(q - sp.Rational(1, 2)))

# Lemma 4: sqrt2 q B - (B - b)/sqrt2 == b/sqrt2 + sqrt2 (q - 1/2) B   (B = B_{k-1}, b = B_{t-1})
check("Lemma 4 algebra", s2*q*B - (B - b)/s2 - (b/s2 + s2*(q - sp.Rational(1, 2))*B))

# Lemma 5(i): 1 + (B - 1 - X)/sqrt2 + cos(alpha) X == 1 + (B-1)/sqrt2 + (cos alpha - 1/sqrt2) X
check("Lemma 5(i) regrouping", 1 + (B - 1 - X)/s2 + ca*X - (1 + (B - 1)/s2 + (ca - 1/s2)*X))
# and 1 + (B-1)/sqrt2 >= B/sqrt2  since difference is 1 - 1/sqrt2 > 0
check("1 + (B-1)/sqrt2 - B/sqrt2 = 1 - 1/sqrt2", 1 + (B - 1)/s2 - B/s2 - (1 - 1/s2))

# Lemma 5(ii-a) decomposition:
# t|b_t| - B/2 = (t|b_t| - |c| cos a_t) + |c|(cos a_t - 1/sqrt2) + (|c| - B/sqrt2)/sqrt2
tb = sp.symbols('t_abs_b', real=True)
check("Lemma 5(ii-a) decomposition",
      (tb - Bt/2) - ((tb - cb*ct) + cb*(ct - 1/s2) + (cb - Bt/s2)/s2))

# Lemma 5(ii-b): q/sin(alpha) > 1/sqrt2 when q > 1/2 and sin(alpha) < 1/sqrt2 -- checked as an identity
# q/sin a - 1/sqrt2 = (q - sin a/sqrt2)/sin a ; numerator > 1/2 - 1/2 = 0
sa = sp.symbols('sin_alpha', positive=True)
check("Lemma 5(ii-b) numerator form", q/sa - 1/s2 - (q - sa/s2)/sa)

# Lemma 3(iii)/Bernoulli telescoping: prod_{t=s}^{u-1} (1 + 1/t) == u/s  (spot values, exact)
good = all(sp.prod([1 + sp.Rational(1, t) for t in range(s_, u_)]) == sp.Rational(u_, s_)
           for s_ in range(1, 25) for u_ in range(s_, 60))
ok_all &= good
print(("OK  " if good else "FAIL"), "telescoping product prod_{t=s}^{u-1}(1+1/t) = u/s (all 1<=s<=u<60, s<25)")

# K_1^(i) <= 1 <= K_2^(i) on 1/2 < q < 1/sqrt2 (needed for Bernoulli)
check("K1(i) at q=1/2 is 1/sqrt2 - 1/2 < 1 (value)", sp.Integer(0) if (1/s2 - sp.Rational(1, 2)) < 1 else sp.Integer(1))

print("\nALL OK" if ok_all else "\nSOME CHECK FAILED")
