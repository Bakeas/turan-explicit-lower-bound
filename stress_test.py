"""
stress_test.py -- test the explicit inequalities of the note on REAL systems.

Theorem 3.1 is a proof by contradiction, so its final step can never be tested
on an actual system. But every intermediate inequality of Appendix A holds for
EVERY system z_1 = 1, ..., z_n with max_j |S_j| <= q, for ANY q in (1/2, 1/sqrt2)
(not only q close to 1/2). Near-optimal systems with max |S_j| ~ 0.69-0.707 < 1/sqrt2
exist (Cheer-Goldston; Bakeas's certificates). So we compute, on such systems, both
sides of every inequality with our explicit constants, and record the smallest
slack. A negative slack would show an error. This is a test, NOT a proof:
it cannot see errors that appear only for q very close to 1/2.
mpmath, 200 digits by default (environment variable DPS); at 60 digits the
polynomial expansion for the n = 300 float tuple loses too many digits.

Usage (systems from https://github.com/Bakeas/turan-power-sum-certificates):
  python stress_test.py small_n_tuples.json results/certs/n_300.json \
      results/certs/n_500.json results/certs/n_1000.json
"""
import json, sys, random
import mpmath as mp

import os
mp.mp.dps = int(os.environ.get("DPS", "200"))
PI = mp.pi
S2 = mp.sqrt(2)
RES = {}      # name -> (min slack, where)


def rec(name, slack, where, rel_to=None):
    s = slack if rel_to in (None, 0) else slack / abs(rel_to)
    if name not in RES or s < RES[name][0]:
        RES[name] = (s, where)


def ang(w1, w2):
    """a(w1, w2) = |arg(w2/w1)|, or pi if w1 w2 = 0."""
    if w1 == 0 or w2 == 0:
        return PI
    return abs(mp.arg(w2 / w1))


def load_small(key):
    d = json.load(open(SMALL))[key]
    z = [mp.mpc(mp.mpf(a), mp.mpf(b)) for a, b in d["z"]]
    n = len(z)
    p = [mp.mpc(1)]
    for r in z[1:]:
        new = p + [mp.mpc(0)]
        for i in range(1, len(new)):
            new[i] = new[i] - r * p[i - 1]
        p = new
    S = [None] + [sum(zz ** j for zz in z) for j in range(1, n + 1)]
    return "small n=%s" % key, n, p, S


def load_cert(path):
    d = json.load(open(path))
    n = d["n"]
    b = [mp.mpc(mp.mpf(int(re[0])) / mp.mpf(2) ** re[1], mp.mpf(int(im[0])) / mp.mpf(2) ** im[1])
         for re, im in d["b"]]
    w = [b[0]] + [b[i] - b[i - 1] for i in range(1, n)] + [-b[n - 1]]
    s = [None] * (n + 1)
    for k in range(1, n + 1):
        acc = -k * w[k]
        for i in range(1, k):
            acc -= w[i] * s[k - i]
        s[k] = acc
    return "cert n=%d" % n, n, b, s


def constants(q, alpha, beta0):
    d = q - mp.mpf(1) / 2
    K1 = 1 / S2 - q
    M = (S2 - 1) / (2 * d)
    kap = M ** (-1 / (1 + q))
    f1 = S2 * d / (mp.cos(alpha) - 1 / S2)
    f2b = q / mp.sin(alpha) - 1 / S2
    f2a = mp.sqrt(q ** 2 - mp.sin(alpha) ** 2 / 2) + min(1, q / mp.sin(alpha)) * (mp.cos(alpha) - 1 / S2) + f2b / S2
    f3 = f2a + f1 / (2 * K1)
    f5p = (2 + q) * (f3 + 3 * d)
    f5pp = (2 + q) * (f3 + 3 * d + (S2 * (1 + q) + mp.mpf(1) / 2) * kap)
    f7 = q ** 2 - mp.mpf(1) / 4
    f8 = (f7 + f5p) / (6 * mp.cos(beta0))
    f8k = (f7 + f5pp) / (6 * mp.cos(beta0))
    return dict(d=d, K1=K1, M=M, kap=kap, f1=f1, f2a=f2a, f2b=f2b, f3=f3, f5p=f5p, f5pp=f5pp,
                f7=f7, f8=f8, f8k=f8k)


def test_system(label, n, b, S):
    q = max(abs(S[j]) for j in range(1, n + 1)) * (1 + mp.mpf(10) ** -40)
    if not (q < 1 / S2):
        print("%s: q = %s >= 1/sqrt2, skipped" % (label, mp.nstr(q, 10)))
        return
    half = mp.mpf(1) / 2
    d = q - half
    c, B, C, D = [], [], [], []
    for t in range(n):
        c.append(b[t] + (c[-1] if t else 0))
        B.append(abs(b[t]) + (B[-1] if t else 0))
        C.append(abs(c[t]) + (C[-1] if t else 0))
        D.append(c[t] + (D[-1] if t else 0))
    alph = [None] + [ang(c[t - 1], b[t]) for t in range(1, n)]
    k = next(t for t in range(1, n) if alph[t] > PI / 4)
    beta = [None] + [ang(D[t - 1], c[t]) for t in range(1, n)]
    print("%s: q = %s, k = %d" % (label, mp.nstr(q, 12), k))
    W = label

    # ---- Lemma 2, (14'), (3) on the whole range ----
    for t in range(1, n):
        rec("L2(i) |c_{t-1}-t b_t| <= q B_{t-1}", q * B[t - 1] - abs(c[t - 1] - t * b[t]), W, B[t - 1])
        rec("(14') |2D_{t-1}-t c_t| <= q C_{t-1}", q * C[t - 1] - abs(2 * D[t - 1] - t * c[t]), W, C[t - 1])
        rec("L2(iii) |c_{t-1}-t b_t| >= |c_{t-1}| sin a_t", abs(c[t - 1] - t * b[t]) - abs(c[t - 1]) * mp.sin(alph[t]), W, B[t - 1])
    for t in range(1, k + 1):
        rec("(star) |c_{t-1}| >= (1-1/s2) + B_{t-1}/s2", abs(c[t - 1]) - (1 - 1 / S2) - B[t - 1] / S2, W, B[t - 1])
    rec("alpha_k < pi/2", PI / 2 - alph[k], W)
    rec("L2(iv) |c_{k-1}| < s2 q B_{k-1}", S2 * q * B[k - 1] - abs(c[k - 1]), W, B[k - 1])
    # ---- Lemma 3, Corollary 1 ----
    Ka, Kb = 1 / S2 - q, 1 + q
    Kc, Kd = 1 / (S2 * (1 + q)), 1 + 1 / (1 / S2 - q)
    for t in range(1, k + 1):
        rec("L3(i) lower Ka B_{t-1} <= t|b_t|", t * abs(b[t]) - Ka * B[t - 1], W, B[t - 1])
        rec("L3(i) upper t|b_t| <= Kb B_{t-1}", Kb * B[t - 1] - t * abs(b[t]), W, B[t - 1])
    for t in range(1, k):
        r = abs(c[t]) / (t * abs(b[t]))
        rec("L3(ii) lower", r - Kc, W)
        rec("L3(ii) upper", Kd - r, W)
    for s_ in range(1, k + 1):
        for u in range(s_, k + 1):
            ratio = B[u - 1] / B[s_ - 1]
            rec("L3(iii) lower (u/s)^Ka", ratio - (mp.mpf(u) / s_) ** Ka, W, ratio)
            rec("L3(iii) upper (u/s)^(1+q)", (mp.mpf(u) / s_) ** Kb - ratio, W, ratio)
    M = (S2 - 1) / (2 * d)
    rec("Cor1 B_{k-1} > M", B[k - 1] - M, W, M)
    rec("Cor1 k > M^(1/(1+q))", k - M ** (1 / (1 + q)), W)
    # ---- Lemma 4 ----
    for t in range(1, k + 1):
        rec("L4 upper |c_{t-1}| <= B_{t-1}/s2 + s2 d B_{k-1}", B[t - 1] / S2 + S2 * d * B[k - 1] - abs(c[t - 1]), W, B[k - 1])
    # ---- Corollary 2 ----
    kap = M ** (-1 / (1 + q))
    rec("1/k < kappa", kap - mp.mpf(1) / k, W)
    rnd = random.Random(1)
    Hs = [list(range(k)), list(range(k // 2 + 1)), list(range(k // 2, k))]
    order = sorted(range(k), key=lambda t: abs(c[t]))
    Hs += [order[:max(1, k // 3)], order[:max(1, k // 10)]]
    Hs += [[t for t in range(k) if rnd.random() < 0.3] or [0] for _ in range(20)]
    for H in Hs:
        y = sum(abs(c[t]) for t in H) / C[k - 1]
        lhs = sum(abs(b[t]) for t in H)
        for eps in (mp.mpf('0.01'), mp.mpf('0.05'), mp.mpf('0.1'), mp.mpf('0.3'), mp.mpf('0.6'), mp.mpf('0.95')):
            bound8 = ((eps + mp.mpf(1) / k) ** Ka + y / (eps * Kc)) * B[k - 1]
            rec("Cor2 (8') with 1/k", bound8 - lhs, W, B[k - 1])
            m4 = min(1, (eps + kap) ** Ka + S2 * (1 + q) * y / eps)
            rec("Cor2 m_4(y) with kappa", m4 * B[k - 1] - lhs, W, B[k - 1])
    # ---- Lemma 5, 6 for several alpha ----
    for alpha in (mp.mpf('0.2'), mp.mpf('0.5'), mp.mpf('0.7'), mp.mpf('0.77'), PI / 4 - mp.mpf('0.001')):
        K = constants(q, alpha, PI / 8)
        Ha = [t for t in range(1, k) if alph[t] <= alpha]
        rec("L5(i) sum_{H_a}|b_t| <= f1 B_{k-1}", K["f1"] * B[k - 1] - sum(abs(b[t]) for t in Ha), W, B[k - 1])
        for t in range(1, k):
            if t in Ha:
                continue
            rec("L5(ii-b) |c_{t-1}| - B_{t-1}/s2 <= f2b B_{t-1}", K["f2b"] * B[t - 1] - (abs(c[t - 1]) - B[t - 1] / S2), W, B[t - 1])
            rec("L5(ii-a) |t|b_t| - B_{t-1}/2| <= f2a B_{t-1}", K["f2a"] * B[t - 1] - abs(t * abs(b[t]) - B[t - 1] / 2), W, B[t - 1])
        for t in range(1, k + 1):
            Pt = sum(j * abs(b[j]) for j in range(1, t + 1))
            At = sum(B[j - 1] for j in range(1, t + 1))
            rec("(13) P_t + A_t = t B_t", -abs(Pt + At - t * B[t]), W, t * B[t])
            if t <= k - 1:
                rec("(11) |P_t - A_t/2| <= t f3 B_{k-1}", t * K["f3"] * B[k - 1] - abs(Pt - At / 2), W, t * B[k - 1])
            rec("(12') 0 <= C_{t-1} - A_t/s2", C[t - 1] - At / S2, W, C[t - 1])
            rec("(12') C_{t-1} - A_t/s2 <= t s2 d B_{k-1}", t * S2 * d * B[k - 1] - (C[t - 1] - At / S2), W, t * B[k - 1])
            rho = t * abs(c[t]) / C[t - 1] - mp.mpf(3) / 2
            X = B[k - 1] / B[t]
            if t <= k - 1:
                rec("L6 |rho_t| <= f5' X_t (t<k)", K["f5p"] * X - abs(rho), W)
                rec("K4 C_{t-1} >= t B_t/(s2(2+q))", C[t - 1] - t * B[t] / (S2 * (2 + q)), W, C[t - 1])
            else:
                rec("L6 |rho_k| <= f5'' X_k (t=k)", K["f5pp"] * X - abs(rho), W)
    # ---- Lemma 7, (17), (18), (19') ----
    bstar = mp.asin(q) / 2
    h = [None] + [abs(D[t - 1]) / C[t - 1] for t in range(1, n)]
    for t in range(1, k + 1):
        rec("L7 beta_t <= pi/8", PI / 8 - beta[t], W)
        rec("L7 (S1) beta_t <= beta* = asin(q)/2", bstar - beta[t], W)
    for t in range(1, k + 2):
        rec("(17) h_t >= cos beta*", h[t] - mp.cos(bstar), W)
        rec("(17) h_t <= 1", 1 - h[t], W)
    F = lambda x: (2 * x * x + 1) / (3 * x)
    for t in range(1, k + 1):
        lam = abs(c[t]) / C[t]
        rec("(18a) h_{t+1} >= h_t + (cos b_t - h_t)|c_t|/C_t", h[t + 1] - (h[t] + (mp.cos(beta[t]) - h[t]) * lam), W)
        if t <= k - 1:
            rec("(18b) (S5) |c_t|/C_t >= 1/(s2 (t+1))", lam - 1 / (S2 * (t + 1)), W, lam)
        rho = t * abs(c[t]) / C[t - 1] - mp.mpf(3) / 2
        f7 = q ** 2 - mp.mpf(1) / 4
        rec("(19') (S2) cos b_t >= F(h_t) - (f7+|rho_t|)/(6 h_t)", mp.cos(beta[t]) - (F(h[t]) - (f7 + abs(rho)) / (6 * h[t])), W)
        for b0 in (PI / 8, bstar):
            K = constants(q, PI / 4 - mp.mpf('0.001'), b0)
            if t <= k - 1:
                rec("(19') with f_8 X_t", mp.cos(beta[t]) - (F(h[t]) - K["f8"] * B[k - 1] / B[t]), W)
            else:
                rec("(19') at t=k with f_8^(k)", mp.cos(beta[t]) - (F(h[t]) - K["f8k"]), W)
                g9 = 1 - mp.cos(b0)          # fallback branch (L8a/L8b fail here)
                rec("L8(ii) fallback cos b_k >= 1 - (g9/3 + f_8^(k))", mp.cos(beta[k]) - (1 - (g9 / 3 + K["f8k"])), W)
    # ---- Section 4 of the paper, with the actual values ----
    hk, Cc, Dd = h[k], C[k - 1], D[k - 1]
    rhok = k * abs(c[k]) / Cc - mp.mpf(3) / 2
    lhs38 = abs(sum(c[j] * S[k - j] for j in range(k)) - Dd / 2)
    e38 = mp.mpf(3) / 2 * (1 - hk) + mp.mpf(3) / 2 * mp.sqrt(2 * (1 - mp.cos(beta[k]))) + abs(rhok)
    rec("(38') actual: |sum c_j S_{k-j} - D/2| <= [3/2(1-h)+3/2 sqrt(2(1-cos b_k))+|rho_k|] C", e38 * Cc - lhs38, W, Cc)
    for b0 in (PI / 8, bstar):
        K = constants(q, PI / 4 - mp.mpf('0.001'), b0)
        g9 = 1 - mp.cos(b0)
        f10 = g9 / 3 + K["f8k"]
        f11 = mp.mpf(3) / 2 * g9 + mp.mpf(3) / 2 * mp.sqrt(2 * f10) + K["f5pp"]
        rec("(38') with f_11 (fallback g_9)", f11 * Cc - lhs38, W, Cc)
    tj = [c[j] / Dd for j in range(k)]
    sabs = sum(abs(x) for x in tj)
    rec("(40') sum|t_j| = 1/h", -abs(sabs - 1 / hk), W)
    s41 = sum(abs(abs(x) - x) for x in tj)
    rec("(41') sum||t_j|-t_j| <= sqrt(2(1-h))/h", mp.sqrt(2 * (1 - hk)) / hk - s41, W)
    z42 = sum(abs(x) * S[k - j] for j, x in enumerate(tj))
    e42 = e38 / hk + q * s41
    rec("(42') |sum|t_j|S - 1/2| <= E38/h + q sum||t|-t|", e42 - abs(z42 - mp.mpf(1) / 2), W)
    f15a = q * sabs - mp.mpf(1) / 2 + e42
    rhs3 = sum(b[j] * S[k - j] for j in range(k))
    rec("(3) at t=k: sum b_j S_{k-j} = c_{k-1} - k b_k", -abs(rhs3 - (c[k - 1] - k * b[k])), W, B[k - 1])
    for gam in (mp.mpf('0.001'), mp.mpf('0.01'), mp.mpf('0.05'), mp.mpf('0.1'), mp.mpf('0.2'), mp.mpf('0.4')):
        cc = mp.mpf(1) / 2 - gam
        Hc = [j for j in range(k) if mp.re(S[k - j]) <= cc]
        Tc = sum(abs(tj[j]) for j in Hc)
        rec("Hc claim (q-c) T_c <= f15 (actual)", f15a - (q - cc) * Tc, W)
        f18 = mp.sqrt(q ** 2 - mp.mpf(1) / 4 + gam)
        for j in range(k):
            if j not in Hc:
                rec("Claim A |S_{k-j} - 1/2| < f18 off H_c", f18 - abs(S[k - j] - mp.mpf(1) / 2), W)
        lhsB = abs(rhs3 - c[k - 1] / 2)
        rhsB = f18 * sum(abs(b[j]) for j in range(k) if j not in Hc) + (q + mp.mpf(1) / 2) * sum(abs(b[j]) for j in Hc)
        rec("Claim B", rhsB - lhsB, W, B[k - 1])
        y = sum(abs(c[j]) for j in Hc) / Cc
        for eps in (mp.mpf('0.05'), mp.mpf('0.3')):
            m4 = min(1, (eps + kap) ** Ka + S2 * (1 + q) * y / eps)
            rec("(43) sum_{H_c}|b_j| <= m_4(y) B", m4 * B[k - 1] - sum(abs(b[j]) for j in Hc), W, B[k - 1])
    e44 = abs(c[k - 1] / 2 - k * b[k])
    rec("(S8) contrapositive: |c_{k-1}/2 - k b_k| > B_{k-1}/4", e44 - B[k - 1] / 4, W, B[k - 1])
    rec("Biro (45)-(47) contrapositive: > (2-s2)/4 B_{k-1}", e44 - (2 - S2) / 4 * B[k - 1], W, B[k - 1])


if __name__ == "__main__":
    SMALL, CERTS = sys.argv[1], sys.argv[2:]
    for key in ("40", "50", "55", "70", "100", "150", "200", "300"):
        test_system(*load_small(key))
    for p in CERTS:
        test_system(*load_cert(p))
    print("\nSmallest slack per inequality (relative where scaled; must be >= 0, > 0 for strict):")
    bad = 0
    for name, (s, where) in RES.items():
        flag = "" if s >= -mp.mpf(10) ** -40 else "   <-- VIOLATED"
        bad += bool(flag)
        print("  %-80s %12s  (%s)%s" % (name, mp.nstr(s, 4), where, flag))
    print("\n%d inequalities tested; %d violated." % (len(RES), bad))
