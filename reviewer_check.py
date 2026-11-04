#!/usr/bin/env python3
"""
reviewer_check.py -- Arb ball-arithmetic checker (python-flint) for the explicit
version of the proof in
  A. Biro, An improved estimate in a power sum problem of Turan,
  Indag. Math. (N.S.) 11 (2000) 343-358.
The method is Biro's; this script only checks the explicit conditions C1-C11
of reviewer_constants.md. Unrefereed.

Parameters (exact, positive):
  eps   = q - 1/2        (q is the constant eta that gets certified)
  a     = pi/4 - alpha
  Lam   = Lambda         (threshold in Lemma 8; replaces f_8^(-1/3))
  phi   = f_9            (Lemma 8)
  gam   = 1/K_6          (Lemma 8)
  theta = auxiliary parameter of Corollary 2
  mu    = 1/2 - c        (Section 4; replaces f_16)
Syntax: name=VALUE, VALUE = p | p/q | 2^(e) | p/q*2^(e)  (e any integer; '**' ok).
On Windows put each argument in double quotes, e.g. "eps=2^(-1000)".
Prints PASS and exits with code 0 only if every condition holds for the whole
ball; otherwise prints FAIL and exits with code 1.
"""
import sys
from flint import arb, ctx

ctx.prec = 256
LOG2 = arb(2).log()

DEFAULTS = {
    "eps":   "2^(-327000000000000000300)",
    "a":     "2^(-218000000000000000160)",
    "Lam":   "2^(109000000000000000000)",
    "phi":   "13/4000000000000000000000",
    "gam":   "1/1000",
    "theta": "1/111000",
    "mu":    "1/716",
}


def parse(text):
    """'p', 'p/q', '2^(e)', 'p/q*2^(e)' -> (p, q, e) meaning (p/q)*2^e."""
    t = text.replace(" ", "").replace("**", "^")
    coef, e = t, 0
    if "^" in t:
        if "*" in t:
            coef, pw = t.rsplit("*", 1)
        else:
            coef, pw = "1", t
        base, ex = pw.split("^", 1)
        if base != "2":
            raise ValueError("only powers of 2 are allowed: " + text)
        e = int(ex.strip("()"))
    if "/" in coef:
        p, q = coef.split("/", 1)
    else:
        p, q = coef, "1"
    p, q = int(p), int(q)
    if p <= 0 or q <= 0:
        raise ValueError("parameters must be positive: " + text)
    return p, q, e


def pow2(e):
    """Exact 2^e as an arb (binary exponentiation; every product is exact)."""
    neg = e < 0
    e = -e if neg else e
    r, b = arb(1), arb(2)
    while e:
        if e & 1:
            r = r * b
        e >>= 1
        if e:
            b = b * b
    return arb(1) / r if neg else r


def value(tup):
    p, q, e = tup
    return (arb(p) / arb(q)) * pow2(e)


def logvalue(tup):
    p, q, e = tup
    return (arb(p) / arb(q)).log() + arb(e) * LOG2


def fmt(x, n=12):
    try:
        return x.str(n, radius=False)
    except Exception:
        return str(x)


def main(argv):
    P = dict(DEFAULTS)
    for arg in argv:
        if "=" not in arg:
            raise ValueError("use name=value, got: " + arg)
        k, v = arg.split("=", 1)
        if k not in P:
            raise ValueError("unknown parameter: " + k)
        P[k] = v
    T = {k: parse(v) for k, v in P.items()}
    eps, a, Lam = value(T["eps"]), value(T["a"]), value(T["Lam"])
    phi, gam = value(T["phi"]), value(T["gam"])
    theta, mu = value(T["theta"]), value(T["mu"])

    one, two = arb(1), arb(2)
    half = one / 2
    r2 = two.sqrt()
    s = one / r2                      # 1/sqrt(2)
    c8 = (two + r2).sqrt() / 2        # cos(pi/8)
    q = half + eps                    # never used inside a difference q - 1/2

    u01 = arb(0.5, 0.5)               # ball containing [0, 1]
    if not (u01.mid() == arb(0.5) and u01.rad() >= arb(0.5)):
        raise RuntimeError("arb(mid, rad) did not produce a ball containing [0,1]")

    # Taylor enclosures, valid for all a >= 0 (notes, Lemma 1 section):
    # sin x in [x - x^3/6, x],  1 - cos x in [x^2/2 - x^4/24, x^2/2]
    a2 = a * a
    sa = a - (a2 * a / 6) * u01                   # sin(a)
    ca1 = a2 / 2 - (a2 * a2 / 24) * u01           # 1 - cos(a)
    b2 = 2 * a
    s2a = b2 - (b2 * b2 * b2 / 6) * u01           # sin(2a)

    # cancellation-free forms, alpha = pi/4 - a  (notes, Lemma 5)
    cosA_s = s * (sa - ca1)                       # cos(alpha) - s
    sinA = s * (1 - ca1 - sa)                     # sin(alpha)
    K1 = (s - half) - eps                         # K_1 = s - q
    K2 = 1 + q                                    # K_2
    kap = s / (1 + q)                             # kappa (Lemma 3 (ii))
    Blow = (r2 - 1) / (2 * eps)                   # Corollary 1
    R = ((r2 - 1) / eps).sqrt() - 1               # k > R
    rho = 1 / R                                   # 1/k < rho

    f1 = r2 * eps / cosA_s                        # Lemma 5 (i)
    dA = (eps + (sa + ca1) / 2) / sinA            # delta_A = q/sin(alpha) - s
    rootA = (eps + eps * eps + s2a / 4).sqrt()    # sqrt(q^2 - sin(alpha)^2/2)
    f2a = cosA_s + rootA + s * dA                 # Lemma 5 (ii)
    E11 = f2a + f1 / (2 * K1) + K2 * rho          # (11)
    E12 = dA + (1 - s) * f1 / K1 + r2 * eps * rho # (12)
    f4 = E11 + 3 * s * E12
    E5 = r2 * eps + (1 + s) * (1 + q) * rho
    f5 = (2 + q) * (r2 * E5 + f4)                 # Lemma 6
    f7 = eps * (1 + eps)                          # q^2 - 1/4

    K7 = two / 3 - gam                            # Lemma 8
    logQlow = logvalue(T["Lam"]) - (2 + q).log()  # log(Lam/(2+q)), valid under C3
    psi = phi / 3 + (f7 + f5) / (6 * c8)          # f_10
    t38 = 3 * phi / 2 + f5
    E38 = (t38 * t38 + 3 * (arb(3) / 2 + f5) * psi).sqrt()     # f_11
    f14 = (E38 + q * (2 * phi).sqrt()) / (1 - phi)
    f15 = eps + q * phi / (1 - phi) + f14
    y = f15 / mu
    m4 = ((theta + rho).log() * K1).exp() + y / (kap * theta)  # f_17 = m_4(y)
    f18 = (mu + f7).sqrt()
    f19 = (q + half) * m4 + f18
    target = (2 - r2) / 4                         # s*(s - 1/2)

    lhs6, rhs6 = f7 + f5 * Lam, 6 * c8 * gam * phi
    lhs7, rhs7 = s * K7 * phi * logQlow, (1 - c8) * (1 + q)

    conds = [
        ("C1  0 < eps <= (sqrt2-1)/4", eps > 0 and eps <= (r2 - 1) / 4, ""),
        ("C2  0 < a <= 1/2, cos(alpha)-s > 0, sin(alpha) > 0",
         a > 0 and a <= half and cosA_s > 0 and sinA > 0, ""),
        ("C3  1 < Lam <= B_low", Lam > 1 and Lam <= Blow, ""),
        ("C4  0 < phi <= 1 - s", phi > 0 and phi <= 1 - s, ""),
        ("C5  0 < gam < 2/3", gam > 0 and gam < two / 3, ""),
        ("C6  f7 + f5*Lam <= 6 cos(pi/8) gam phi", lhs6 <= rhs6,
         "log2(rhs/lhs) = " + fmt((rhs6 / lhs6).log() / LOG2, 6)),
        ("C7  s K7 phi log(Lam/(2+q)) > (1-cos(pi/8))(1+q)", lhs7 > rhs7,
         "lhs/rhs = " + fmt(lhs7 / rhs7, 8)),
        ("C8  0 < theta <= 1", theta > 0 and theta <= 1, ""),
        ("C9  mu > 0", mu > 0, ""),
        ("C10 0 < K1, 2 K1 <= 1, R >= 1", K1 > 0 and 2 * K1 <= 1 and R >= 1, ""),
        ("C11 f19 < (2-sqrt2)/4", f19 < target,
         "f19 = " + fmt(f19, 10) + ", slack = " + fmt(target - f19, 6)),
    ]

    print("python-flint Arb check, prec =", ctx.prec)
    for k in ("eps", "a", "Lam", "phi", "gam", "theta", "mu"):
        print("  %-5s = %s" % (k, P[k]))
    print("  log2(f5) =", fmt(f5.log() / LOG2, 15))
    print("  f15 =", fmt(f15, 10), "  y =", fmt(y, 10), "  m4 =", fmt(m4, 10))
    ok = True
    for name, res, info in conds:
        res = bool(res)
        ok = ok and res
        print(("  ok   " if res else "  FAIL ") + name
              + (("   [" + info + "]") if info else ""))
    if ok:
        print("PASS")
        print("certified (given the derivations in reviewer_constants.md):")
        print("  max_{1<=j<=n} |S_j| > 1/2 + eps for every n and all z_1 = 1, z_2, ..., z_n,")
        print("  with eps =", P["eps"])
        return 0
    print("FAIL")
    return 1


if __name__ == "__main__":
    try:
        code = main(sys.argv[1:])
    except Exception as exc:
        print("FAIL (error: %r)" % (exc,))
        code = 1
    sys.exit(code)