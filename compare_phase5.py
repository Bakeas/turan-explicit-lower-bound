# Phase 5 numeric comparison: our chain (check_eta.compute)
# versus the second checker's chain (reviewer_check.py), at matched parameters.
# Only a consistency test, NOT a proof. The proofs are compared in Section 6 of the note.
import sys
from flint import arb, ctx
import check_eta as ours

ctx.prec = 256
LOG2 = arb(2).log()


def rev_values(args):
    src = open("reviewer_check.py").read()
    marker = '    print("python-flint Arb check, prec =", ctx.prec)'
    assert marker in src
    src = src.replace(marker, "    return locals()\n" + marker)
    g = {"__name__": "rev"}
    exec(compile(src, "reviewer_check.py", "exec"), g)
    return g["main"](args)


def our_values(**kw):
    P = {
        "delta": ours.Param(kw["delta"]), "x": ours.Param(kw["x"]),
        "theta": ours.Param(kw["theta"], allow_one_minus=True),
        "mu": ours.Param(kw["mu"]), "eps": ours.Param(kw["eps"]),
        "gamma": ours.Param(kw["gamma"]), "beta0": kw["beta0"],
    }
    V, C, finb = ours.compute(P)
    d = {n: b for n, b in V}
    ok = all(st is True for _, _, st in C)
    return d, C, finb, ok


def s(x, n=10):
    return x.str(n, radius=False)


# 1. sanity: our chain at the point of Theorem 1.1
d, C, finb, ok = our_values(delta="2^-32176926823186696", x="2^-21451284548791130",
                            beta0="star", theta="1-2^-40", mu="2^-30",
                            eps="12/100000", gamma="408/100000")
print("[1] Theorem 1.1 point, our chain: verdict", "PASS" if ok else "FAIL",
      "| f_19 =", s(d["f_19"], 20))

# 2. the second checker's point (run 2) in both chains
N_eps, N_a, N_L = 318600000000000000400, 212400000000000000200, 106200000000000000000
R = rev_values(["eps=2^(-%d)" % N_eps, "a=2^(-%d)" % N_a, "Lam=2^(%d)" % N_L,
                "phi=33/10000000000000000000000", "gam=1/100000",
                "theta=1/110000", "mu=1/714"])
# map: our delta=eps, x=a, mu(L8)=gam, eps(Cor2)=theta, gamma(c)=mu,
# beta0=pi/8; our theta chosen so that f_8^(-theta) = Lam, i.e.
# theta = log2(Lam)/log2(1/f_8).
d0, _, _, _ = our_values(delta="2^-%d" % N_eps, x="2^-%d" % N_a, beta0="pi8",
                         theta="1/2", mu="1/100000", eps="1/110000", gamma="1/714")
L8 = -(d0["f_8"].log() / LOG2)          # log2(1/f_8), ball
c = L8 - N_L                             # log2(1/(f_8 Lam)) > 0
print("[2] our log2(1/f_8) =", s(L8, 25), "; minus log2(Lam) =", s(c, 10))
# theta = 1 - c/L8; use a rational slightly larger c' >= c (theta a bit smaller,
# so our f_8^(-theta) <= Lam, the same side as the second checker's C6).
cnum = int(c.upper().floor().unique_fmpz()) + 1 if hasattr(c, "upper") else None
den = int(L8.lower().floor().unique_fmpz())
theta_txt = "1-%d/%d" % (cnum, den)
d, C, finb, ok = our_values(delta="2^-%d" % N_eps, x="2^-%d" % N_a, beta0="pi8",
                            theta=theta_txt, mu="1/100000", eps="1/110000",
                            gamma="1/714")
print("    our theta =", theta_txt)
rows = [
    ("f_5 (reviewer) | f_5'' (ours)", R["f5"], d["f_5'' (t = k)"]),
    ("phi (reviewer, chosen) | g_9 (ours, computed)", R["phi"], d["g_9 = min(f_9, 1 - cos beta_0)"]),
    ("f_10: psi | f_10", R["psi"], d["f_10"]),
    ("f_11: E38 | f_11", R["E38"], d["f_11"]),
    ("f_14", R["f14"], d["f_14"]),
    ("f_15", R["f15"], d["f_15"]),
    ("y", R["y"], d["y = f_15/gamma"]),
    ("m4(y) | f_17", R["m4"], d["f_17 = min(1, ...)"]),
    ("f_18", R["f18"], d["f_18 = sqrt(f_7 + gamma)"]),
    ("f_19", R["f19"], d["f_19"]),
]
print("    %-48s %-20s %-20s %s" % ("quantity", "reviewer", "ours", "ours/reviewer"))
for name, a, b in rows:
    print("    %-48s %-20s %-20s %s" % (name, s(a), s(b), s(b / a, 8)))
print("    Biro's threshold (2-sqrt2)/4 =", s((2 - arb(2).sqrt()) / 4, 12))
print("    our chain at this point: conditions",
      {lab: st for lab, _, st in C}, "| FIN-B:", finb)
