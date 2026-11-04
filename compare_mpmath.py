#!/usr/bin/env python3
"""
compare_mpmath.py -- cross-check for check_eta.py (NOT a proof).

An independent float prototype of explicit_constants.md section 2.27, written
with mpmath at very high precision, using the DIRECT formulas (q, alpha,
1/sqrt2 - q, q/sin(alpha) - 1/sqrt2, ...), not the cancellation-free forms.
At moderate points (q - 1/2 >= 1e-300) enough digits are carried that the
direct forms are accurate.  It compares every value with the midpoint
computed by check_eta.compute() and reports the largest relative difference.
"""
import sys
from fractions import Fraction
import mpmath as mp

import check_eta as ce


def proto(delta, x, beta0, theta, mu, eps, gamma, cf=False):
    """cf=False: direct formulas (needs many digits);
    cf=True: cancellation-free forms of section 2.27 (any size)."""
    q = mp.mpf(1) / 2 + delta
    s2 = mp.sqrt(2)
    if not cf:
        al = mp.pi / 4 - x
        d = q - mp.mpf(1) / 2
        K1 = 1 / s2 - q
        sa, cam = mp.sin(al), mp.cos(al) - 1 / s2
        f2b = q / sa - 1 / s2
        r2a = q ** 2 - sa ** 2 / 2
        f7 = q ** 2 - mp.mpf(1) / 4
    else:
        d = delta
        K1 = (1 / s2 - mp.mpf(1) / 2) - d
        sa = (mp.cos(x) - mp.sin(x)) / s2
        cam = (mp.sin(x) - 2 * mp.sin(x / 2) ** 2) / s2
        f2b = (d + mp.sin(x / 2) ** 2 + mp.sin(x) / 2) / sa
        r2a = d + d ** 2 + mp.sin(2 * x) / 4
        f7 = d * (1 + d)
    M = (s2 - 1) / (2 * d)
    kap = mp.exp(-mp.log(M) / (1 + q))
    f1 = s2 * d / cam
    f2a = mp.sqrt(r2a) + min(1, q / sa) * cam + f2b / s2
    f3 = f2a + f1 / (2 * K1)
    f5p = (2 + q) * (f3 + 3 * d)
    f5pp = (2 + q) * (f3 + 3 * d + (s2 * (1 + q) + mp.mpf(1) / 2) * kap)
    b0 = mp.pi / 8 if beta0 == "pi8" else mp.asin(q) / 2
    f8 = (f7 + f5p) / (6 * mp.cos(b0))
    f8k = (f7 + f5pp) / (6 * mp.cos(b0))
    K11 = s2 * (1 + q) * (1 - mp.cos(b0)) / (mp.mpf(2) / 3 - mu)
    m5 = M / (2 + q)
    K12 = 1 / (2 + q)
    f9 = max(f8 ** (1 - theta) / mu, K11 / mp.log(m5),
             K11 / mp.log(K12 * f8 ** (-theta)))
    g9 = min(f9, 1 - mp.cos(b0))
    f10 = g9 / 3 + f8k
    f11 = mp.mpf(3) / 2 * g9 + mp.mpf(3) / 2 * mp.sqrt(2 * f10) + f5pp
    f12 = g9 / (1 - g9)
    f13 = mp.sqrt(2 * g9) / (1 - g9)
    f14 = f11 / (1 - g9) + q * f13
    f15 = d + q * f12 + f14
    y = f15 / gamma
    f17 = min(1, (eps + kap) ** K1 + s2 * (1 + q) * y / eps)
    f18 = mp.sqrt(f7 + gamma)
    f19 = f18 + (q + mp.mpf(1) / 2) * f17
    return {"K_1 = 1/sqrt2 - q": K1, "kappa = M^(-1/(1+q))": kap, "f_1": f1,
            "f_2a": f2a, "f_3": f3, "f_5' (t <= k-1)": f5p,
            "f_5'' (t = k)": f5pp, "f_8": f8, "f_8^(k)": f8k, "K_11": K11,
            "log m_5": mp.log(m5), "f_9": f9, "f_10": f10, "f_11": f11,
            "f_12": f12, "f_13": f13, "f_14": f14, "f_15": f15,
            "f_17 = min(1, ...)": f17, "f_18 = sqrt(f_7 + gamma)": f18,
            "f_19": f19}


def q2mp(p):
    v = mp.mpf(p.num) / p.den * mp.mpf(2) ** p.exp
    return (1 - v) if p.one_minus else v


def run(args, cf=False):
    a = dict(args)
    P = {k: ce.Param(a[k], allow_one_minus=(k == "theta"))
         for k in ("delta", "x", "theta", "mu", "eps", "gamma")}
    P["beta0"] = a["beta0"]
    ce.ctx.prec = 1024 if not cf else 256
    V, C, _ = ce.compute(P)
    arbv = dict(V)
    pv = proto(*(q2mp(P[k]) for k in ("delta", "x")), a["beta0"],
               *(q2mp(P[k]) for k in ("theta", "mu", "eps", "gamma")), cf=cf)
    worst = 0
    for k, v in pv.items():
        m = mp.mpf(arbv[k].mid().str(250, radius=False).replace(" ", ""))
        rel = abs(m - v) / abs(v) if v != 0 else abs(m)
        worst = max(worst, rel)
    return worst, arbv["f_19"].str(8, radius=False), [s for _, _, s in C]


if __name__ == "__main__":
    mp.mp.dps = 1200
    moderate = [
        dict(delta="1/1000000", x="1/1000", beta0="pi8", theta="49/100",
             mu="3/10", eps="1/100", gamma="1/100"),
        dict(delta="2^-40", x="2^-20", beta0="star", theta="44/100",
             mu="19/100", eps="1/1000", gamma="1/50"),
        dict(delta="3*2^-200", x="5*2^-100", beta0="star", theta="3/4",
             mu="1/40", eps="1/10000", gamma="1/200"),
        dict(delta="1/" + str(10 ** 400), x="2^-660",
             beta0="pi8", theta="1-2^-10", mu="1/200", eps="1/100000",
             gamma="1/300"),
        # uncapped branch of f_17 (m_4 < 1): large gamma, eps
        dict(delta="2^-40", x="2^-20", beta0="pi8", theta="1/2",
             mu="1/3", eps="1/2", gamma="100"),
    ]
    huge = [
        dict(delta="2^-2308312065422341452", x="2^-1442695040888963407",
             beta0="star", theta="1-2^-52", mu="1/1000", eps="28/1000000",
             gamma="25/10000"),
        dict(delta="2^-23083120654223415", x="2^-14426950408889634",
             beta0="star", theta="1-2^-48", mu="1/1000", eps="12/100000",
             gamma="45/10000"),
    ]
    bad = False

    def show(tag, p, worst, f19):
        print("%-6s %-44s max rel diff %.1e  f_19 = %s" % (
            tag, ",".join("%s=%s" % (k, v if len(v) < 10 else v[:8] + "..")
                          for k, v in p.items() if k in ("delta", "x", "beta0", "gamma")),
            worst, f19))

    for p in moderate:
        for cf in (False, True):
            if cf:
                mp.mp.dps = 60
            else:
                mp.mp.dps = 1200
            worst, f19, conds = run(p, cf=cf)
            show("cf" if cf else "direct", p, worst, f19)
            if worst > (mp.mpf(10) ** -50 if cf else mp.mpf(10) ** -100):
                bad = True
    mp.mp.dps = 60
    for p in huge:
        worst, f19, conds = run(p, cf=True)
        show("huge", p, worst, f19)
        if worst > mp.mpf(10) ** -40:
            bad = True
    print("COMPARISON", "FAILED" if bad else "OK")
    sys.exit(1 if bad else 0)
