#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_eta.py  --  Phase 3 checker for an explicit lower bound in Turan's
pure power sum problem (Erdos problem #519, constant C42).

WHAT IT DOES
  Given the parameters (q, alpha, beta_0, theta, mu, eps, gamma) of the
  explicit Theorem in explicit_constants.md, section 2.27, it computes the
  constants f_1, ..., f_19 (steps 1-15 of section 2.27) in Arb ball
  arithmetic (python-flint) and checks the conditions
      (P1)-(P4), (L8a), (L8b), (FIN)
  for the whole ball, i.e. rigorously.  (FIN-B), Biro's original final
  threshold, is printed for information only.
  It prints PASS and exits with code 0 only if every condition holds.
  Otherwise it prints FAIL and exits with code 1 (also for bad input or
  any error).

  If it prints PASS, then (assuming the written derivations in
  explicit_constants.md are correct; they are UNREFEREED)
      max_{1<=j<=n} |z_1^j + ... + z_n^j| > q
  for every n and all complex z_1 = 1, z_2, ..., z_n, with q = 1/2 + delta.

CREDIT
  The method and all ideas are A. Biro's: "An improved estimate in a power
  sum problem of Turan", Indag. Math. (N.S.) 11 (2000), 343-358.  We only
  make his constants explicit (and add a few small sharpenings inside his
  argument, listed in explicit_constants.md section 5).
  Prepared by V. Bakeas with Claude (AI).  Unrefereed working tool.

PARAMETERS (all exact)
  q     = 1/2 + delta,   delta = a/b * 2^e          (--delta)
  alpha = pi/4 - x,      x     = a/b * 2^e          (--x)
          (alpha is only an auxiliary parameter; it need not be rational)
  beta_0 = pi/8  or  (1/2) arcsin q                 (--beta0 pi8 | star)
  theta  = a/b * 2^e   or   1 - a/b * 2^e           (--theta)
  mu, eps, gamma = a/b * 2^e                        (--mu --eps --gamma)
  Number format: "a/b*2^e", "a*2^e", "a/b", "a", or "2^e"
  (a, b positive integers, e any integer, e.g. "3*2^-45", "1/1000",
  "2^-40000000000000000").  For theta also "1-2^-45" etc.
  Huge exponents (like 2^-4e16) are fine: Arb exponents are arbitrary
  size, and only cancellation-free formulas are used (section 2.27).

EXAMPLE
  python check_eta.py --delta "2^-20" --x "2^-10" --beta0 pi8 \
      --theta 49/100 --mu 3/10 --eps 1/100 --gamma 1/100

REQUIREMENTS
  python-flint (tested with 0.9.0).
"""

import argparse
import json
import re
import sys
from fractions import Fraction

try:
    from flint import arb, ctx, fmpq, fmpz
except ImportError:  # pragma: no cover
    print("ERROR: python-flint is not installed (pip install python-flint).")
    print("FAIL")
    sys.exit(1)


# --------------------------------------------------------------------------
# Input parsing (exact)
# --------------------------------------------------------------------------

_NUM_RE = re.compile(
    r"^\s*(?:(?P<num>\d+)(?:\s*/\s*(?P<den>\d+))?)?"
    r"\s*(?:\*?\s*2\s*\^\s*\(?\s*(?P<exp>[-+]?\d+)\s*\)?)?\s*$"
)


class Param:
    """A positive exact number a/b * 2^e, optionally in the form 1 - a/b*2^e."""

    def __init__(self, text, allow_one_minus=False):
        self.text = text.strip()
        s = self.text.replace(" ", "")
        self.one_minus = False
        if allow_one_minus and s.startswith("1-"):
            self.one_minus = True
            s = s[2:]
        m = _NUM_RE.match(s)
        if not m or (m.group("num") is None and m.group("exp") is None):
            raise ValueError("cannot parse number %r (use a/b*2^e)" % text)
        num = int(m.group("num")) if m.group("num") is not None else 1
        den = int(m.group("den")) if m.group("den") is not None else 1
        exp = int(m.group("exp")) if m.group("exp") is not None else 0
        if num <= 0 or den <= 0:
            raise ValueError("number %r must be positive" % text)
        self.num, self.den, self.exp = num, den, exp

    def core(self):
        """Ball containing a/b * 2^e (exact if b is a power of 2)."""
        return arb(fmpq(self.num, self.den)) * arb((fmpz(1), fmpz(self.exp)))

    def value(self):
        """Ball containing the parameter value."""
        c = self.core()
        return (1 - c) if self.one_minus else c

    def complement(self):
        """Ball containing 1 - value (exact-form when given as 1 - c)."""
        c = self.core()
        return c if self.one_minus else (1 - c)

    def core_cmp(self, r):
        """Exact sign of (a/b*2^e - r) for a positive Fraction r."""
        A = self.num * r.denominator
        B = r.numerator * self.den
        e = self.exp
        if e >= 0:
            if e > B.bit_length():
                return 1
            lhs, rhs = A << e, B
        else:
            if -e > A.bit_length():
                return -1
            lhs, rhs = A, B << (-e)
        return (lhs > rhs) - (lhs < rhs)

    def exact_str(self):
        core = "%d" % self.num
        if self.den != 1:
            core += "/%d" % self.den
        if self.exp != 0:
            core += "*2^(%d)" % self.exp
        return ("1 - " + core) if self.one_minus else core


# --------------------------------------------------------------------------
# Output helpers
# --------------------------------------------------------------------------

def fmt(v, digits=12):
    """Readable midpoint of a ball; huge exponents are fine."""
    try:
        if not v.is_finite():
            return "not finite"
        return v.str(digits, radius=False)
    except Exception:  # pragma: no cover
        return str(v)


def acc(v):
    b = v.rel_accuracy_bits()
    return "exact" if b > 10**9 else "%d bits" % b


# --------------------------------------------------------------------------
# The computation of section 2.27
# --------------------------------------------------------------------------

def compute(P):
    """Steps 1-15 of explicit_constants.md section 2.27, in Arb balls.

    Returns (values, conditions).  values: ordered list of (name, ball).
    conditions: ordered list of (label, statement, status) where status is
    True (holds for the whole ball), False (fails for the whole ball) or
    None (undecided at this precision: counts as failure).
    """
    V = []      # (name, ball)
    C = []      # (label, statement, status)

    def val(name, b):
        V.append((name, b))
        return b

    def cond(label, text, true_if, false_if):
        st = True if true_if else (False if false_if else None)
        C.append((label, text, st))
        return st

    one = arb(1)
    half = arb(1) / 2
    s2 = arb(2).sqrt()
    inv_s2 = 1 / s2
    pi = arb.pi()

    # ----- parameters ---------------------------------------------------
    delta = val("delta = q - 1/2", P["delta"].value())
    q = val("q = 1/2 + delta (eta)", half + delta)
    x = val("x = pi/4 - alpha", P["x"].value())
    theta = val("theta", P["theta"].value())
    one_m_theta = val("1 - theta", P["theta"].complement())
    mu = val("mu", P["mu"].value())
    eps = val("eps", P["eps"].value())
    gamma = val("gamma", P["gamma"].value())

    # ----- (P1)-(P4) --------------------------------------------------------
    gap = inv_s2 - half            # 1/sqrt2 - 1/2
    cond("(P1)", "1/2 < q < 1/sqrt2  [i.e. 0 < delta < 1/sqrt2 - 1/2]",
         (delta > 0) and (delta < gap), (delta <= 0) or (delta >= gap))
    cond("(P2)", "0 < alpha < pi/4  [i.e. 0 < x < pi/4]",
         (x > 0) and (x < pi / 4), (x <= 0) or (x >= pi / 4))
    C.append(("(P3)", "beta_0 = %s" % (
        "pi/8" if P["beta0"] == "pi8" else "(1/2) arcsin q"), True))
    two3 = arb(2) / 3
    # (P4) is decided exactly (rational arithmetic).  All inputs a/b*2^e are
    # > 0 by parsing, so: 0 < theta < 1  <=>  core(theta) < 1 (both forms);
    # mu < 2/3; eps < 1; gamma > 0 is automatic.
    p4 = (P["theta"].core_cmp(Fraction(1)) < 0
          and P["mu"].core_cmp(Fraction(2, 3)) < 0
          and P["eps"].core_cmp(Fraction(1)) < 0)
    C.append(("(P4)", "0 < theta < 1, 0 < mu < 2/3, 0 < eps < 1, gamma > 0"
              "  [exact]", p4))

    # ----- step 1 -----------------------------------------------------------
    K1 = val("K_1 = 1/sqrt2 - q", gap - delta)
    M = val("M = (sqrt2 - 1)/(2 delta)", (s2 - 1) / (2 * delta))
    logM = M.log()
    kappa = val("kappa = M^(-1/(1+q))", (-logM / (1 + q)).exp())

    # ----- trigonometric quantities, cancellation-free (G1-G4) -----------
    sx, cx = x.sin(), x.cos()
    sh2 = (x / 2).sin() ** 2                      # sin^2(x/2) = (1 - cos x)/2
    sin_a = val("sin alpha = (cos x - sin x)/sqrt2", (cx - sx) * inv_s2)
    cos_a_m = val("cos alpha - 1/sqrt2 = (sin x - 2 sin^2(x/2))/sqrt2",
                  (sx - 2 * sh2) * inv_s2)
    rad2a = val("q^2 - sin^2(alpha)/2 = delta + delta^2 + sin(2x)/4",
                delta + delta ** 2 + (2 * x).sin() / 4)

    # ----- step 2 -----------------------------------------------------------
    f1 = val("f_1", s2 * delta / cos_a_m)

    # ----- step 3 -----------------------------------------------------------
    f2b = val("f_2b = (delta + sin^2(x/2) + sin(x)/2)/sin alpha",
              (delta + sh2 + sx / 2) / sin_a)
    min1 = (q / sin_a).min(one)
    f2a = val("f_2a", rad2a.sqrt() + min1 * cos_a_m + f2b * inv_s2)

    # ----- step 4 -----------------------------------------------------------
    f3 = val("f_3", f2a + f1 / (2 * K1))

    # ----- step 5 -----------------------------------------------------------
    f5p = val("f_5' (t <= k-1)", (2 + q) * (f3 + 3 * delta))
    f5pp = val("f_5'' (t = k)",
               (2 + q) * (f3 + 3 * delta + (s2 * (1 + q) + half) * kappa))

    # ----- step 6 -----------------------------------------------------------
    f7 = val("f_7 = q^2 - 1/4 = delta(1 + delta)", delta * (1 + delta))
    if P["beta0"] == "pi8":
        beta0 = pi / 8
    else:
        beta0 = q.asin() / 2
    beta0 = val("beta_0", beta0)
    cb0 = val("cos beta_0", beta0.cos())
    omc = val("1 - cos beta_0 = 2 sin^2(beta_0/2)", 2 * (beta0 / 2).sin() ** 2)
    f8 = val("f_8", (f7 + f5p) / (6 * cb0))
    f8k = val("f_8^(k)", (f7 + f5pp) / (6 * cb0))
    logf8 = f8.log()

    # ----- step 7 -----------------------------------------------------------
    K11 = val("K_11", s2 * (1 + q) * omc / (two3 - mu))
    m5 = val("m_5 = M/(2 + q)", M / (2 + q))
    K12 = val("K_12 = 1/(2 + q)", 1 / (2 + q))

    # ----- (L8a), (L8b) -------------------------------------------------------
    log_m5 = val("log m_5", logM - (2 + q).log())
    L8a = cond("(L8a)", "m_5(q) > 1", log_m5 > 0, log_m5 <= 0)
    # log(K_12 f_8^(-theta)) = log K_12 - theta log f_8
    log_K12f8 = val("log(K_12 f_8^(-theta))", K12.log() - theta * logf8)
    lhs_L8b = val("(2 + q) f_8^theta", (2 + q) * (theta * logf8).exp())
    L8b = cond("(L8b)", "(2 + q) f_8^theta < 1",
               (lhs_L8b < 1) and (log_K12f8 > 0),
               (lhs_L8b >= 1) or (log_K12f8 <= 0))

    # ----- step 8 -----------------------------------------------------------
    if L8a is True and L8b is True:
        t1 = val("f_9 term 1: f_8^(1-theta)/mu",
                 (one_m_theta * logf8).exp() / mu)
        t2 = val("f_9 term 2: K_11/log m_5", K11 / log_m5)
        t3 = val("f_9 term 3: K_11/log(K_12 f_8^(-theta))", K11 / log_K12f8)
        f9 = val("f_9", t1.max(t2).max(t3))
        g9 = val("g_9 = min(f_9, 1 - cos beta_0)", f9.min(omc))
    else:
        # section 2.20: g_9 := 1 - cos beta_0 if (L8a) or (L8b) fails.
        # (Here only for information; the verdict is FAIL anyway.)
        g9 = val("g_9 = 1 - cos beta_0  [(L8a)/(L8b) not verified]", omc)

    # ----- step 9-12 ----------------------------------------------------------
    f10 = val("f_10", g9 / 3 + f8k)
    f11 = val("f_11", 3 * g9 / 2 + 3 * (2 * f10).sqrt() / 2 + f5pp)
    f12 = val("f_12", g9 / (1 - g9))
    f13 = val("f_13", (2 * g9).sqrt() / (1 - g9))
    f14 = val("f_14", f11 / (1 - g9) + q * f13)
    f15 = val("f_15", delta + q * f12 + f14)

    # ----- step 13 ------------------------------------------------------------
    y = val("y = f_15/gamma", f15 / gamma)
    m4a = val("(eps + kappa)^K_1", (K1 * (eps + kappa).log()).exp())
    m4b = val("sqrt2 (1 + q) y/eps", s2 * (1 + q) * y / eps)
    f17 = val("f_17 = min(1, ...)", (m4a + m4b).min(one))

    # ----- step 14 ------------------------------------------------------------
    f18 = val("f_18 = sqrt(f_7 + gamma)", (f7 + gamma).sqrt())
    f19 = val("f_19", f18 + (q + half) * f17)

    # ----- consistency checks (proved in the text from (P1), (P2)) ------------
    cond("(S1)", "K_1 > 0, sin alpha > 0, cos alpha - 1/sqrt2 > 0",
         (K1 > 0) and (sin_a > 0) and (cos_a_m > 0),
         (K1 <= 0) or (sin_a <= 0) or (cos_a_m <= 0))
    cond("(S2)", "q^2 - sin^2(alpha)/2 > 0, f_2b > 0",
         (rad2a > 0) and (f2b > 0), (rad2a <= 0) or (f2b <= 0))
    cond("(S3)", "cos beta_0 > 1/sqrt2", cb0 > inv_s2, cb0 <= inv_s2)
    cond("(S4)", "0 < g_9 < 1", (g9 > 0) and (g9 < 1), (g9 <= 0) or (g9 >= 1))

    # ----- step 15 ------------------------------------------------------------
    quarter = arb(1) / 4
    cond("(FIN)", "f_19 <= 1/4", f19 <= quarter, f19 > quarter)
    finb = (2 - s2) / 4
    finb_st = True if (f19 <= finb) else (False if (f19 > finb) else None)
    return V, C, finb_st


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

class _Parser(argparse.ArgumentParser):
    def error(self, message):
        self.print_usage(sys.stderr)
        print("ERROR: %s" % message)
        print("FAIL")
        sys.exit(1)


def main(argv=None):
    ap = _Parser(description="Rigorous checker for the explicit lower bound "
                 "(explicit_constants.md section 2.27).")
    ap.add_argument("--delta", required=True, help="q - 1/2, as a/b*2^e")
    ap.add_argument("--x", required=True, help="pi/4 - alpha, as a/b*2^e")
    ap.add_argument("--beta0", required=True, choices=["pi8", "star"],
                    help="pi8: beta_0 = pi/8;  star: beta_0 = (1/2) arcsin q")
    ap.add_argument("--theta", required=True, help="a/b*2^e or 1-a/b*2^e")
    ap.add_argument("--mu", required=True)
    ap.add_argument("--eps", required=True)
    ap.add_argument("--gamma", required=True)
    ap.add_argument("--prec", type=int, default=256,
                    help="working precision in bits (default 256)")
    ap.add_argument("--json", help="also write all values to this JSON file")
    ap.add_argument("--quiet", action="store_true",
                    help="print only the conditions and the verdict")
    a = ap.parse_args(argv)

    try:
        P = {
            "delta": Param(a.delta),
            "x": Param(a.x),
            "theta": Param(a.theta, allow_one_minus=True),
            "mu": Param(a.mu),
            "eps": Param(a.eps),
            "gamma": Param(a.gamma),
            "beta0": a.beta0,
        }
        if a.prec < 64:
            raise ValueError("--prec must be at least 64")
        ctx.prec = a.prec
        V, C, finb_st = compute(P)
    except Exception as e:
        print("ERROR: %s: %s" % (type(e).__name__, e))
        print("FAIL")
        return 1

    print("check_eta.py -- explicit form of Biro's lower bound [Bir00b]")
    print("(conditions of explicit_constants.md section 2.27; unrefereed)")
    print("precision: %d bits" % a.prec)
    print()
    print("Parameters (exact):")
    print("  q = eta = 1/2 + %s" % P["delta"].exact_str())
    print("  alpha   = pi/4 - %s" % P["x"].exact_str())
    print("  beta_0  = %s" % ("pi/8" if a.beta0 == "pi8" else "(1/2) arcsin q"))
    for k in ("theta", "mu", "eps", "gamma"):
        print("  %-7s = %s" % (k, P[k].exact_str()))
    print()
    if not a.quiet:
        print("Values (ball midpoints, 12 digits; accuracy = relative bits):")
        for name, b in V:
            print("  %-58s %-26s [%s]" % (name, fmt(b), acc(b)))
        print()

    print("Conditions (each must hold for the whole ball):")
    ok = True
    for label, text, st in C:
        word = "ok" if st is True else ("FAILS" if st is False
                                        else "UNDECIDED (raise --prec)")
        print("  %-6s %-58s %s" % (label, text, word))
        if st is not True:
            ok = False
    w = ("holds" if finb_st is True else
         ("does not hold" if finb_st is False else "undecided"))
    print("  info:  (FIN-B) f_19 <= (2 - sqrt2)/4 (Biro's final step): %s" % w)
    print()

    if a.json:
        try:
            data = {
                "parameters": {k: (P[k].exact_str() if k != "beta0" else P[k])
                               for k in P},
                "prec": a.prec,
                "values": {n: {"mid": fmt(b, 20), "ball": b.str(20)}
                           for n, b in V},
                "conditions": {l: st for l, _, st in C},
                "FIN-B": finb_st,
                "verdict": "PASS" if ok else "FAIL",
            }
            with open(a.json, "w", encoding="utf-8") as fh:
                json.dump(data, fh, indent=1)
        except Exception as e:
            print("ERROR writing JSON: %s" % e)
            ok = False

    if ok:
        print("PASS: all conditions hold. If the derivations in "
              "explicit_constants.md are correct (unrefereed),")
        print("then R_n > eta = 1/2 + %s for every n." % P["delta"].exact_str())
        return 0
    print("FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())
