# reviewer_constants.md

> **Note added by V. Bakeas (2026-11-04).** These are the notes of the separate AI instance that wrote `reviewer_check.py`, working only from Biró's papers [Bir00b] and [B1], without access to our derivations or our checker. They are reproduced unchanged below this box, except for this note. Its notation differs from the note's: its ε, a, P_t, C_t, D_t are the note's δ, x, B_t, D_t, C_t. Its erratum [E1] in §14 is **not** an erratum: on the printed page, Lemma 1(iii) of [Bir00b] reads "≥" (the ">" came from a low-resolution scan). Its [E2] is the one real misprint. Its "Run 2" in §16 is the point of Theorem 1.2 of the note.

Explicit constants for A. Biró, "An improved estimate in a power sum problem of Turán",
Indag. Math. (N.S.) 11 (2000) 343–358 [Bir00b], with A. Biró, Acta Math. Hungar. 65
(1994) 209–216 [B1].

Status: UNREFEREED. The method is Biró's; we only make his constants explicit.
Steps not in [Bir00b] are tagged [S1], [S2], ... (list in §13). Errata in [Bir00b]
are tagged [E1], ... (§14). Every item is marked [derived] or [open].

reviewer_check.py (python-flint 0.9.0, Arb, prec 256) printed PASS for
  eps = 2^(-327000000000000000300), a = 2^(-218000000000000000160),
  Lam = 2^(109000000000000000000), phi = 13/4000000000000000000000,
  gam = 1/1000, theta = 1/111000, mu = 1/716.
If the derivations below are correct: max_{1<=j<=n} |S_j| > eta := 1/2 + eps.

## 0. Notation, standing assumption

- q = 1/2 + ε, s = 1/√2, α = π/4 − a, c8 = cos(π/8) = √(2+√2)/2.
- Assumption (1): |S_j| ≤ q for 1 ≤ j ≤ n. We derive a contradiction from (1)
  and C1–C11 (§12). n = 1 is trivial (|S_1| = 1 > q), so n ≥ 2.
- b_0 = c_0 = 1, c_t = 1 + b_1 + … + b_t.
- P_t = Σ_{j=0}^{t} |b_j| = 1 + |b_1| + … + |b_t|; non-decreasing; |c_t| ≤ P_t.
- C_{t-1} = Σ_{j=0}^{t-1} c_j,  D_{t-1} = Σ_{j=0}^{t-1} |c_j|.
- α_t, β_t, H_α as in [Bir00b] (α_t, β_t := π when undefined).
- Checker names: K1 = s − q, K2 = 1 + q, κ = s/(1+q) (kap),
  B_low = (√2−1)/(2ε), R = √((√2−1)/ε) − 1, ρ = 1/R.
- Bernoulli: for 0 ≤ K ≤ 1 and x ≥ 0, (1+x)^K ≤ 1 + Kx (concavity of x ↦ x^K).

## 1. Source check: [B1] → [Bir00b] (3), (4)   [derived]

[B1] p. 210: (z−z_2)…(z−z_n) = z^{n−1} + b_1 z^{n−2} + … + b_{n−1},
T_j = Σ_{t=2}^{n} z_t^j. Newton's identities (b_j = (−1)^j e_j, degree n−1):
T_k + b_1T_{k−1} + … + b_{k−1}T_1 + k b_k = 0 (1 ≤ k ≤ n−1) and
T_n + b_1T_{n−1} + … + b_{n−1}T_1 = 0 (k = degree + 1, no k b_k term).
With T_j = S_j − 1 (z_1 = 1) this gives [B1](1), (2).
[Bir00b](3) = [B1](1) with k renamed t, same range 1 ≤ t ≤ n−1; [Bir00b](4) = [B1](2).
Signs, index ranges and the b_j agree. [Bir00b] uses (3) only for t ≤ k ≤ n−1
and (4) only for the existence of k. (14) follows from (3) as stated in [Bir00b]
(t = 1: S_1 = 2 − c_1 = 1 − b_1; differences: Σ (c_i − c_{i−1}) S_{t−i} = Σ b_i S_{t−i}).

## 2. Lemma 1   [derived]

(i), (ii), (iv), (v) are correct as stated. (iv): |w1+w2| ≥ Re((w1+w2) conj(w1)/|w1|)
= |w1| + |w2| cos φ.
[E1] (iii) should read |w1 − w2| ≥ |w1|·|sin φ| (equality when |w2| = |w1| cos φ).
It follows from (i). Only the ≥ form is used later.
Iterated (iv): |c_m| ≥ |c_{m−1}| + |b_m| cos α_m for all m ≥ 1. This is trivial if
c_{m−1} = 0 or b_m = 0, since then α_m = π and cos π = −1.
Taylor enclosures used by the checker (a ≥ 0):
  a − a³/6 ≤ sin a ≤ a,   a²/2 − a⁴/24 ≤ 1 − cos a ≤ a²/2,
  cos α − s = s(sin a − (1 − cos a)),   sin α = s(1 − (1 − cos a) − sin a),
  q² − sin²α/2 = ε + ε² + sin(2a)/4   (sin²α = (1 − sin 2a)/2).

## 3. Lemma 2 and the index k   [derived]

(i) For 1 ≤ t ≤ n−1: |c_{t−1} − t b_t| ≤ q P_{t−1}. Proof: (3), (1), triangle inequality.
(ii) |c_{t−1}| ≥ 1 + Σ_{j<t} |b_j| cos α_j (iterated Lemma 1(iv)). With (5):
  (L2*)  |c_{t−1}| ≥ 1 + s(P_{t−1} − 1) ≥ s P_{t−1}   for 1 ≤ t ≤ k.
Existence of k: if α_t ≤ π/4 for all t ≤ n−1, then (L2*) at t = n gives
|c_{n−1}| ≥ sP_{n−1}, while (4) and (1) give |c_{n−1}| ≤ qP_{n−1}. This contradicts
q < s. Here q_0 := (1+√2)/4, and C1 gives q ≤ q_0 < s.
α_k < π/2 and b_k ≠ 0: c_{k−1} ≠ 0 by (L2*). Suppose b_k = 0 or α_k ≥ π/2. Then
|c_{k−1} − k b_k| ≥ |c_{k−1}| (Lemma 1(iv) with w2 = −k b_k, angle π − α_k ≤ π/2),
so sP_{k−1} ≤ qP_{k−1}, a contradiction.
(iii) |c_{t−1} − t b_t| ≥ |c_{t−1}| sin α_t (by [E1]; trivial if α_t = π).
(iv) At t = k: |c_{k−1}| sin α_k ≤ qP_{k−1} and sin α_k > s, so
  |c_{k−1}| < √2 q P_{k−1}   (strict).

## 4. Lemma 3, Corollary 1   [derived]

Lemma 3(i), K1 = s − q, K2 = 1 + q: for 1 ≤ t ≤ k, K1 P_{t−1} ≤ t|b_t| ≤ K2 P_{t−1}.
Proof: t|b_t| ≥ |c_{t−1}| − |c_{t−1} − t b_t| ≥ (s − q)P_{t−1} by (L2*) and 2(i).
Also t|b_t| ≤ |c_{t−1}| + qP_{t−1} ≤ (1+q)P_{t−1}. (K1 > 0 is C10.)
  (P)  (1 + K1/t) P_{t−1} ≤ P_t ≤ (1 + K2/t) P_{t−1} ≤ (2+q) P_{t−1},  1 ≤ t ≤ k.
Lemma 3(ii): for 1 ≤ t ≤ k−1, κ ≤ |c_t|/(t|b_t|) ≤ (2+q)/K1, with κ = s/(1+q).
Proof: |c_t| ≥ sP_t ≥ sP_{t−1} ((L2*) at t+1 ≤ k), t|b_t| ≤ (1+q)P_{t−1}.
Upper bound: |c_t| ≤ P_t ≤ (2+q)P_{t−1}, t|b_t| ≥ K1 P_{t−1}. Only κ is used later.
m2(k) = k(k+1)/2: K2 < 2, so P_{k−1} ≤ Π_{t=1}^{k−1} (t+2)/t = k(k+1)/2.
m1(k) = k^{K1} (from (P) and Bernoulli, K1 ≤ 1). Not used.
Corollary 1 (m3): (L2*) at t = k and 2(iv) give 1 + s(P_{k−1} − 1) < √2 q P_{k−1},
i.e. 1 − s < √2 ε P_{k−1}, so
  P_{k−1} > (1 − s)/(√2 ε) = B_low.
Then B_low < k(k+1)/2 < (k+1)²/2, so (k+1)² > (√2−1)/ε, i.e. k > R.
With R ≥ 1 (C10): k ≥ 2 and 1/k < ρ.

## 5. Corollary 2 (m4)   [derived]

Statement. Let θ ∈ (0, 1] (C8), 0 < K1 ≤ 1 (C10), y ≥ 0. Assume
H ⊆ {0, …, k−1} and Σ_{t∈H} |c_t| ≤ y D_{k−1}. Then
  Σ_{t∈H} |b_t| ≤ m4(y) P_{k−1},   m4(y) = (θ + ρ)^{K1} + y/(κθ).
Proof. Let m = min(⌊θk⌋, k−1).
(a) t ∈ H, t ≤ m: these terms sum to at most P_m (b_0 = 1 is included). By (P) and
Bernoulli: P_{k−1}/P_m ≥ Π_{t=m+1}^{k−1} (1 + K1/t) ≥ Π ((t+1)/t)^{K1} = (k/(m+1))^{K1}
(empty product if m = k−1). Since (m+1)/k ≤ θ + 1/k < θ + ρ, we get
P_m ≤ (θ+ρ)^{K1} P_{k−1}.
(b) t ∈ H, t > m: then 1 ≤ t ≤ k−1 and t > ⌊θk⌋, so t > θk. By Lemma 3(ii),
|b_t| ≤ |c_t|/(κt) < |c_t|/(κθk). Sum ≤ y D_{k−1}/(κθk) ≤ y P_{k−1}/(κθ),
since D_{k−1} ≤ k P_{k−1}.
Relation to [Bir00b]: ε → θ, K3 → K1, the y/(εK1) term → y/(κθ). "k arbitrarily
large" is replaced by the explicit ρ term.

## 6. Lemma 4 and |c_t|   [derived]

For 1 ≤ t ≤ k:  0 ≤ |c_{t−1}| − s P_{t−1} ≤ √2 ε P_{k−1}   (O-constant √2).
Proof: the lower bound is (L2*). For the upper bound, iterate Lemma 1(iv) from t to
k−1 with α_j ≤ π/4: |c_{k−1}| ≥ |c_{t−1}| + s(P_{k−1} − P_{t−1}). Combine with 2(iv).
Consequence (E5): for 1 ≤ t ≤ k, ||c_t| − sP_t| ≤ E5 P_{k−1}, where
E5 = √2ε + (1+s)(1+q)ρ.
For t ≤ k−1 this is Lemma 4 at t+1. For t = k: ||c_k| − |c_{k−1}|| ≤ |b_k| and
sP_k = sP_{k−1} + s|b_k|, so |c_k| − sP_k lies in
[(|c_{k−1}| − sP_{k−1}) − (1+s)|b_k|, (|c_{k−1}| − sP_{k−1}) + (1−s)|b_k|].
Finally k|b_k| ≤ K2 P_{k−1} and 1/k < ρ.

## 7. Lemma 5   [derived]

(i) f1 = √2 ε/(cos α − s), needs cos α > s (C2).
Proof: |c_{k−1}| ≥ 1 + Σ_{j=1}^{k−1} |b_j| cos α_j ≥ sP_{k−1} + (cos α − s) Σ_{H_α} |b_j|.
Compare with 2(iv).
(ii) Let 1 ≤ t ≤ k−1, t ∉ H_α. Then α < α_t ≤ π/4, so sin α_t > sin α > 0 and
0 ≤ cos α_t − s < cos α − s. Write c = c_{t−1}, P = P_{t−1}.
(9): from 2(i) and Lemma 1(i) (w1 = c, w2 = t b_t),
  q²P² ≥ |c|² sin²α_t + (t|b_t| − |c| cos α_t)².
Second formula (checker dA):
  0 ≤ |c| − sP ≤ δ_A P,   δ_A = q/sin α − s = (ε + (sin a + 1 − cos a)/2)/sin α.
(Proof: |c| ≤ qP/sin α_t by (9). For the closed form use s·sin α = (1 − (1−cos a) − sin a)/2.)
First formula (checker f2a):
  |t|b_t| − P/2| ≤ f2 P,   f2 = (cos α − s) + √(q² − sin²α/2) + s δ_A.
Proof: by (9) and |c| ≥ sP, |t|b_t| − |c| cos α_t| ≤ P√(q² − sin²α_t/2) ≤ P√(q² − sin²α/2).
Also |c| cos α_t − P/2 = |c|(cos α_t − s) + s(|c| − sP), which lies in
[0, (cos α − s)P + sδ_A P].
[Bir00b] centres at s|c| instead of P/2; the content is the same.

## 8. Lemma 6   [derived]

Statement: for 1 ≤ t ≤ k,  | t|c_t| / D_{t−1} − 3/2 | ≤ f5 P_{k−1}/P_t, where
  E11 = f2 + f1/(2K1) + K2 ρ,
  E12 = δ_A + (1−s) f1/K1 + √2 ε ρ,
  f4  = E11 + 3s·E12,
  f5  = (2+q)(√2·E5 + f4).
Needs 2K1 ≤ 1 (C10).
Proof. Fix 1 ≤ t ≤ k. Let X = Σ_{j=1}^{t} j|b_j|, Y = Σ_{j=1}^{t} P_{j−1},
Z = Σ_{j=1}^{t} |c_{j−1}| = D_{t−1}.
(13) X + Y = tP_t, by induction: going from t to t+1 both sides grow by
(t+1)|b_{t+1}| + P_t.
(11) |X − Y/2| ≤ t E11 P_{k−1}. Sum |j|b_j| − P_{j−1}/2| over j ≤ t:
  - j ∉ H_α, j ≤ k−1: each term ≤ f2 P_{k−1} (Lemma 5(ii)).
  - j ∈ H_α: both quantities are ≥ 0, so the term is ≤ max(j|b_j|, P_{j−1}/2)
    ≤ j|b_j|/(2K1) (Lemma 3(i), 2K1 ≤ 1). And Σ_{j∈H_α, j≤t} j|b_j| ≤ t f1 P_{k−1}.
  - j = k (only if t = k): term ≤ max(k|b_k|, P_{k−1}/2) ≤ K2 P_{k−1} < t K2 ρ P_{k−1}.
(12) |Z − sY| ≤ t E12 P_{k−1}. Each term |c_{j−1}| − sP_{j−1} is ≥ 0 by (L2*):
  - j ∉ H_α, j ≤ k−1: ≤ δ_A P_{k−1}.
  - j ∈ H_α: ≤ (1−s) P_{j−1} ≤ (1−s) j|b_j|/K1.
  - j = k: ≤ √2 ε P_{k−1} < t √2 ε ρ P_{k−1} (Lemma 4).
Combine: (3/2)Y = tP_t − (X − Y/2). Since 3s·s = 3/2,
  3sZ = (3/2)Y + 3s(Z − sY) = tP_t + e1,   |e1| ≤ t f4 P_{k−1}.
By §6, |c_t| = sP_t + e2 with |e2| ≤ E5 P_{k−1}.
Lower bound on 3sZ: X ≤ K2 Y (Lemma 3(i)), so tP_t ≤ (2+q)Y. Also Z ≥ sY by (L2*),
so 3sZ ≥ (3/2)Y ≥ (3/2) tP_t/(2+q). (Paper's K4 = s/(2+q).)
Then t|c_t|/Z − 3/2 = (3s·t·e2 − (3/2) e1)/(3sZ), so
  |…| ≤ (2+q)(√2|e2| + |e1|/t)/P_t ≤ f5 P_{k−1}/P_t.

## 9. Lemma 7 (β_t ≤ π/8)   [derived]

Needs q < s (C1). For 1 ≤ t ≤ k, (14) and (1) give |2C_{t−1} − t c_t| ≤ q D_{t−1}.
Let t ≤ k be the least index with β_t > β := π/8. Iterating Lemma 1(iv) with
β_j ≤ β for j < t gives (15): |C_{t−1}| ≥ D_{t−1} cos β > 0.
If β_t ≤ π/2: by Lemma 1(iii) [E1], |2C − t c_t| ≥ 2|C| sin β_t ≥ 2|C| sin β.
If β_t > π/2 (this includes c_t = 0): by Lemma 1(iv) with w2 = −t c_t,
|2C − t c_t| ≥ 2|C| ≥ 2|C| sin β.
So q D_{t−1} ≥ 2 sin β cos β D_{t−1} = s D_{t−1}, a contradiction.

## 10. Lemma 8   [derived; S1–S3 are not in Bir00b, see §13]

10.1 Definitions. h_t = |C_{t−1}|/D_{t−1} for 1 ≤ t ≤ k.
(17) c8 ≤ h_t ≤ 1: iterate Lemma 1(iv) with β_j ≤ π/8 (Lemma 7); upper bound by
the triangle inequality.
Step inequality [E2]: for 1 ≤ t ≤ k−1,
  h_{t+1} ≥ (1 − λ_t) h_t + λ_t cos β_t,   λ_t = |c_t|/D_t ∈ (0,1].
Proof: |C_t| ≥ |C_{t−1}| + |c_t| cos β_t (Lemma 1(iv); C_{t−1} ≠ 0 by (17),
|c_t| ≥ sP_t > 0); divide by D_t = D_{t−1} + |c_t|.
Also λ_t ≥ s/(t+1), because |c_t| ≥ sP_t (L2*) and D_t ≤ Σ_{j≤t} P_j ≤ (t+1)P_t.

10.2 Explicit (19) [S1]. Let 1 ≤ t ≤ k, h = h_t, X_t = P_{k−1}/P_t, δ = f5 X_t.
Assume f7 + δ ≤ 3/2. Then
  cos β_t ≥ F(h) − (f7 + δ)/(6h),   F(h) = (2h²+1)/(3h),   f7 = q² − 1/4 = ε(1+ε).
Proof. Put w1 = 2C_{t−1}/D_{t−1} and w2 = t c_t/D_{t−1}. Then |w1| = 2h and
u := |w2| = 3/2 + d with |d| ≤ δ (Lemma 6). Also |w1 − w2| ≤ q by (14) and (1).
Since u ≥ 3/2 − δ ≥ f7 > 0, we have c_t ≠ 0.
Lemma 1(ii): q² ≥ 4h² + u² − 4hu cos β_t, so cos β_t ≥ g := (4h² + u² − q²)/(4hu).
Expanding:
  12hu(g − F(h)) = 3(4h² + u² − q²) − 4u(2h²+1) = −3f7 + d(5 − 8h²) + 3d².
Claim: this is ≥ −(3 + 2d)(f7 + δ) = −2u(f7 + δ), i.e.
  d(5 − 8h²) + 3d² + 3δ + 2d(f7 + δ) ≥ 0.
- d ≥ 0: 5 − 8h² ≥ −3 (h ≤ 1), so d(5 − 8h²) ≥ −3δ; all other terms are ≥ 0.
- d < 0: 8h² ≥ 8c8² = 4 + 2√2 > 5, so d(5 − 8h²) > 0; and
  3δ + 2d(f7 + δ) ≥ δ(3 − 2(f7 + δ)) ≥ 0.
Divide by 12hu > 0.
Where it is used:
- t = k: X_k ≤ 1, so δ ≤ f5.
- t0 ≤ t ≤ k−1: X_t ≤ Q ≤ Λ, so by (17) and C6, (f7 + δ)/(6h) ≤ (f7 + f5Λ)/(6c8) ≤ γφ.
The hypothesis f7 + δ ≤ 3/2 holds because f7 + f5Λ ≤ 6c8γφ < 4(1 − s) < 3/2
(by C4–C6) and f5 ≤ f5Λ (since Λ > 1).

10.3 Facts about F. On (0,1], F(h) − h = (1−h²)/(3h) is decreasing.
F'(h) = (2h²−1)/(3h²) ≥ 0 on [s,1]. For 0 < φ ≤ 1−s [S3]:
  F(1−φ) = 1 − φ(1−2φ)/(3(1−φ)) ≥ 1 − φ/3,
  F(1−φ) − (1−φ) = φ(2−φ)/(3(1−φ)) ≥ 2φ/3.

10.4 Choice of t0 [S2]. Let t0 be the least t ∈ {1,…,k−1} with P_{k−1}/P_t ≤ Λ.
It exists: k ≥ 2, and t = k−1 qualifies since 1 < Λ (C3).
Q := P_{k−1}/P_{t0} ≤ Λ.
- If t0 ≥ 2: P_{k−1}/P_{t0−1} > Λ by minimality.
- If t0 = 1: P_{k−1}/P_0 = P_{k−1} > B_low ≥ Λ (Corollary 1 and C3).
In both cases P_{t0} ≤ (2+q)P_{t0−1} by (P), so Q > Λ/(2+q).

10.5 Proof of (i).
Case A: h_t ≥ 1−φ for some t ∈ [t0, k−1]. Then h_t ∈ [1−φ, 1] ⊂ [s,1] (C4), so
cos β_t ≥ F(1−φ) − γφ ≥ 1 − φ/3 − γφ ≥ 1 − φ (γ ≤ 2/3, C5).
By the step inequality h_{t+1} ≥ 1−φ. Repeating up to k gives h_k ≥ 1−φ.
Case B: h_t < 1−φ for all t ∈ [t0, k−1]. Then
cos β_t − h_t ≥ F(1−φ) − (1−φ) − γφ ≥ K7φ > 0, with K7 = 2/3 − γ.
With λ_t ≥ s/(t+1):
  1 − c8 ≥ h_k − h_{t0} ≥ s K7 φ Σ_{t=t0}^{k−1} 1/(t+1) ≥ s K7 φ Σ_{m=t0+1}^{k−1} 1/m.
By (P), log Q = Σ_{m=t0+1}^{k−1} log(P_m/P_{m−1}) ≤ K2 Σ_{m=t0+1}^{k−1} 1/m. Hence
  (1 − c8)(1+q) ≥ s K7 φ log Q > s K7 φ log(Λ/(2+q)),
which contradicts C7. So Case A holds:
  (i)  h_k ≥ 1 − φ        (f9 = φ).
(ii) Apply 10.2 at t = k with h_k ∈ [1−φ, 1]:
  cos β_k ≥ F(1−φ) − (f7 + f5)/(6c8) ≥ 1 − ψ,   ψ = φ/3 + (f7 + f5)/(6c8)   (f10 = ψ).

## 11. Section 4, (38)–(47)   [derived]

(38) Let D = D_{k−1}. By (14) at t = k, Σ_{j=0}^{k−1} c_j S_{k−j} = 2C_{k−1} − k c_k. So
  |Σ c_j S_{k−j} − C_{k−1}/2| = |(3/2)C_{k−1} − k c_k| ≤ f11 D,
  f11 = E38 = √((3φ/2 + f5)² + 3(3/2 + f5)ψ).
Proof: put v = C/D and w = k c_k/D. Then |v| = h_k ∈ [1−φ, 1] and |w| = 3/2 + d,
|d| ≤ f5 (Lemma 6, t = k). The angle between v and w is β_k. Lemma 1(ii) gives
|(3/2)v − w|² = ((3/2)h_k − |w|)² + 3h_k|w|(1 − cos β_k).
Here |(3/2)(h_k − 1) − d| ≤ 3φ/2 + f5, and 3h_k|w|(1 − cos β_k) ≤ 3(3/2 + f5)ψ.
(39), (40): t_j = c_j/C_{k−1}, Σ t_j = 1, Σ|t_j| = 1/h_k ≤ 1/(1−φ). So f12 = φ/(1−φ).
(41): Lemma 1(v) and Cauchy–Schwarz give Σ ||t_j| − t_j| ≤ √2 √(Σ|t_j| − 1) √(Σ|t_j|)
  ≤ √(2φ)/(1−φ) = f13.
(42): divide (38) by C_{k−1} (the error becomes f11/h_k), then replace t_j by |t_j| (error ≤ q f13):
  |Σ |t_j| S_{k−j} − 1/2| ≤ f14 = (E38 + q√(2φ))/(1−φ).
Take real parts: Σ |t_j| Re S_{k−j} ≥ 1/2 − f14.
Let c = 1/2 − μ (μ > 0, C9; [Bir00b] allows any c < 1/2 and takes μ = √f15) and
H_c = {0 ≤ j ≤ k−1: Re S_{k−j} ≤ c}. Then
Σ|t_j| Re S_{k−j} ≤ qΣ|t_j| − (q − c)Σ_{H_c}|t_j|, and q − c = ε + μ ≥ μ, so
  Σ_{H_c} |t_j| ≤ y := f15/μ,   f15 = ε + qφ/(1−φ) + f14.
(43): Σ_{H_c}|c_j| ≤ y|C_{k−1}| ≤ y D_{k−1}. Corollary 2 with H = H_c gives
  Σ_{H_c}|b_j| ≤ m4(y) P_{k−1}   (f17 = m4(y)).
For j ∉ H_c: Re S > 1/2 − μ and |S| ≤ q, so
|S − 1/2|² = |S|² − Re S + 1/4 < f7 + μ. Hence f18 = √(μ + f7).
For j ∈ H_c: |S − 1/2| ≤ q + 1/2. Therefore
  |Σ_{j=0}^{k−1} b_j S_{k−j} − c_{k−1}/2| ≤ f19 P_{k−1},   f19 = (q + 1/2) m4(y) + f18.
(44): by (3) at t = k ≤ n−1, Σ b_j S_{k−j} = c_{k−1} − k b_k, so
|c_{k−1}/2 − k b_k| ≤ f19 P_{k−1}.
(45)–(47): by 2(iii), sin α_k > s and c_{k−1} ≠ 0, so |c_{k−1} − k b_k| > s|c_{k−1}|. Hence
  (s − 1/2)|c_{k−1}| < |c_{k−1} − k b_k| − |c_{k−1}/2| ≤ f19 P_{k−1},
and with (46) = (L2*), s(s − 1/2) < f19 (f20 = f19/(s − 1/2)).
This contradicts C11. So (1) is false, i.e. max_{1≤j≤n} |S_j| > q.  ∎

## 12. Conditions used by the final contradiction (= checker C1–C11)

C1  0 < ε ≤ (√2−1)/4. Gives q ≤ q0 = (1+√2)/4 < s.
    Used in §3 (k exists, α_k < π/2), §9, and for K1 > 0.
C2  0 < a ≤ 1/2, cos α − s > 0, sin α > 0. Gives 0 < α < π/4. Used in §7.
C3  1 < Λ ≤ B_low. Used in §10.4 and for f5 ≤ f5Λ.
C4  0 < φ ≤ 1 − s. Used in §10 ([1−φ, 1] ⊂ [s,1]) and §11.
C5  0 < γ < 2/3. Used in §10 (Case A, K7 > 0).
C6  f7 + f5Λ ≤ 6 c8 γ φ. Used in §10.2.
C7  s K7 φ log(Λ/(2+q)) > (1 − c8)(1+q). Used in §10.5, Case B.
C8  0 < θ ≤ 1. Used in §5.
C9  μ > 0. Used in §11.
C10 K1 > 0, 2K1 ≤ 1, R ≥ 1. Used in §4, §5 (K1 ≤ 1), §8 (2K1 ≤ 1), ρ = 1/R, k ≥ 2.
C11 f19 < (2−√2)/4 = s(s − 1/2). Used in §11.
Formula chain (same order as in the checker): K1, K2, κ, B_low, R, ρ; f1, δ_A, f2;
E11, E12, f4, E5, f5; f7, K7, ψ; E38, f14, f15, y, m4(y), f18, f19.

12.2 Map of [Bir00b] constants
K1 = s − q (Lemma 3(i)). In Lemma 3(ii) and Corollary 2 the lower constant is κ = s/(1+q).
K2 = 1 + q (Lemma 3(i)). The upper constant of Lemma 3(ii) is (2+q)/K1 (unused).
K3 = K1.  K4 = s/(2+q).  K5 = s/2 (we use s/(t+1)).  K6 = 1/γ.  K7 = 2/3 − γ.
K8 = s·K7 (with 1/(t+1)).  K9 = 1 − c8.  K10 = 1/(1+q).
K11 = (1 − c8)(1+q)/(s K7).  K12 = 1/(2+q).
K13, K14, K15, m5: not needed. The parameters φ, Λ and C3, C7 take their role [S2].
m1(k) = k^{K1}, m2(k) = k(k+1)/2, m3(q) = min(R, B_low), m4(y) = (θ+ρ)^{K1} + y/(κθ).
f1 §7; f2 → (δ_A, f2); f3 → (E11, E12); f4, f5 §8; f6 → (f7 + δ)/(6h) [S1];
f7 = ε(1+ε); f8 → f5·Λ; f9 = φ; f10 = ψ; f11 = E38; f12 = φ/(1−φ);
f13 = √(2φ)/(1−φ); f14; f15; f16 → μ; f17 = m4(y); f18 = √(μ + f7); f19;
f20 = f19/(s − 1/2).
O-constant in Lemma 4: √2. "k arbitrarily large" → the ρ terms.
"y small enough" → not needed (m4 holds for all y ≥ 0).
[Bir00b] has no δ_i; δ_A and δ are our names.

## 13. Steps NOT in [Bir00b] (sharpenings; need expert check)

[S1] §10.2: the error in (19) is (f7 + f5 X_t)/(6h_t), linear in X_t = P_{k−1}/P_t,
     via the exact identity −3f7 + d(5 − 8h²) + 3d². [Bir00b] has O(f8 X_t²).
[S2] §10.4: t0 is defined by a free threshold Λ instead of f8^{−1/3} ((29)–(37)).
     Case (29) cannot occur, because P_{k−1} > B_low ≥ Λ (C3).
[S3] §10.3: we use the exact value of F(1−φ) instead of (20)'s H + (1−H)/3.
     This allows γ < 2/3 in Case A and gives ψ = φ/3 + … in (ii).
Explicit versions of steps that ARE in [Bir00b] (not sharpenings): the ρ terms
(Lemma 6 at t = k, Corollary 2), the Bernoulli product bound, and the free μ in §11.

## 14. Errata in [Bir00b] (read from an OCR'd PDF; confirm on the printed page)

[E1] Lemma 1(iii): should be |w1 − w2| ≥ |w1|·|sin φ|, not ">". Harmless.
[E2] (18): the left side should be h_{t+1}, not h_{t−1}. The preceding display bounds h_{t+1}.

## 15. Open items

- No gap found in [Bir00b].
- O1: expert check of [S1]–[S3]. Full proofs are given above, but they are new.
- O2: confirm E1, E2 on the printed page.
- O3: an independent line-by-line check of notes ↔ reviewer_check.py.

## 16. Optimization and eta

Float model (FLOAT, guidance only, derived by hand):
- F1: f5 ≈ 2.5(√(a/2) + 4.83 ε/a). Balancing gives a ∝ ε^{2/3} and f5 ≍ ε^{1/3},
  so C6 gives log2(1/ε) ≈ 3 log2 Λ + O(1).
- F2: f15 ≈ 1.932√φ and f19 ≈ θ^{K1} + 1.932√φ/(μκθ) + √μ. The optimum is
  θ^{K1} = A/(1 + 3K1) with A = (2−√2)/4, and √μ = 2(A − θ^{K1})/3.
  This gives θ ≈ 9.1e-6, μ ≈ 1/714, φ_max ≈ 3.37e-21.
- F3: C7 gives log2 Λ ≈ 0.3494/φ.
  Floor for this chain of estimates: log2(1/ε) ≈ 3.11e20.
Bottleneck: θ^{K1} with K1 ≈ 0.207 (Corollary 2) forces θ ~ 1e-5, so y ~ 1e-7.
Then φ ~ y² by (41), log Λ ~ 1/φ by Lemma 8, and log(1/ε) ~ 3 log Λ by Lemma 5.

Run 1 (PASS, 2026-11-04, user's machine, python-flint 0.9.0;
self-test arb(0,1) > 0 → False, arb(2,1) > 0 → True):
  eps = 2^(-327000000000000000300), a = 2^(-218000000000000000160),
  Lam = 2^(109000000000000000000), phi = 13/4000000000000000000000,
  gam = 1/1000, theta = 1/111000, mu = 1/716.
  C6: log2(rhs/lhs) = 3.62; C7: lhs/rhs = 1.0122; C11: slack 3.34e-4.
  Result: eta = 1/2 + 2^(-327000000000000000300) ≈ 1/2 + 10^(-9.84e19).
  UNREFEREED: it holds only if §0–§11 are correct.
Run 2 (PASS, 2026-11-04, user's machine, python-flint 0.9.0):
  eps = 2^(-318600000000000000400), a = 2^(-212400000000000000200),
  Lam = 2^(106200000000000000000), phi = 33/10000000000000000000000,
  gam = 1/100000, theta = 1/110000, mu = 1/714.
  C6: log2(rhs/lhs) = 17.0011; C7: lhs/rhs = 1.0029019;
  C11: f19 = 0.1462564428, slack 1.90e-4.
  Best result: eta = 1/2 + 2^(-318600000000000000400) ≈ 1/2 + 10^(-9.59e19).
  UNREFEREED: holds only if §0–§11 are correct (see §13 S1–S3, §15).
