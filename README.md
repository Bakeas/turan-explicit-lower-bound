> **Note (4 Nov 2026).** After the first version was published, I learned that R. Pinho had already proposed, on 2 Nov 2026, a much better explicit value from Biró's method, $C_{42}\ge 1/2+2^{-4096}$* (unrefereed; [teorth/optimizationproblems PR #310](https://github.com/teorth/optimizationproblems/pull/310)). His argument replaces part of Biró's Lemma 8 by a contraction, which gives a much faster rate. The bound in this repository is far weaker: it is the value that Biró's argument gives as it is written. Version 1.0.1 of the note cites his work (Section 4.4).

# An explicit lower bound in Turán's pure power sum problem

**Claim.** For every $n\ge1$ and all complex numbers $z_1=1,z_2,\dots,z_n$,

$$\max_{1\le k\le n}\Bigl|\sum_{j=1}^n z_j^k\Bigr| > \frac12+2^{-32176926823186696}.$$

Hence $R_n>\frac12+2^{-32176926823186696}>\frac12+10^{-10^{16}}$ for every $n$, where $R_n=\min_{z_1=1}\max_{1\le k\le n}\lvert z_1^k+\dots+z_n^k\rvert$. This is Turán's pure power sum problem ([Erdős problem #519](https://www.erdosproblems.com/519); constant [C42](https://teorth.github.io/optimizationproblems/constants/42a.html) of the Optimization Constants repository).

A more conservative value, whose proof uses Biró's own Lemma 7 ($\beta_t\le\pi/8$) and his own final step (45)–(47): $R_n>\frac12+2^{-318600000000000000400}>\frac12+10^{-10^{20}}$.

**Credit.** The method is entirely A. Biró's: *An improved estimate in a power sum problem of Turán*, Indag. Math. (N.S.) **11** (2000), 343–358. Biró proved $R_n>q$ for an effectively computable constant $q>1/2$, but did not compute $q$. We follow his proof step by step, make every constant explicit, and add ten small sharpenings S1–S10 inside his argument; all of them are listed, with proofs, in Section 4 of the note.

**Status.** Unrefereed. No professional mathematician has checked the argument. Made with substantial help from an AI model (see [AI use](#ai-use)). Corrections are very welcome (GitHub issue or email).

**Versions.**
- 1.0.1 (4 Nov 2026): the note cites R. Pinho's earlier and much better explicit value (introduction and Section 4.4). Mathematics and checkers unchanged.
- 1.0.0 (4 Nov 2026): first version, [doi:10.5281/zenodo.23961978](https://doi.org/10.5281/zenodo.23961978).

**Author:** Vasileios Bakeas, Department of Mechanical Engineering, University of Thessaly, Volos, Greece (bakeasb@gmail.com)

## The note

`Bakeas_explicit_lower_bound_Turan.pdf` (LaTeX source: `Bakeas_explicit_lower_bound_Turan.tex`).

- Theorem 3.1: the explicit form of Biró's theorem. If seven parameters satisfy the conditions (P1)–(P4), (L8a), (L8b) and (FIN), then $\max_j\lvert S_j\rvert>q$. All quantities are given by 14 explicit formulas.
- Section 4: every step that is not in Biró's paper (one misprint, the changes S1–S10, bookkeeping).
- Section 5: the two parameter points and the rigorous numerical verification.
- Appendix A: the complete proof, line by line, in the order of Biró's paper.

## Verify

```bash
pip install python-flint==0.9.0

# Theorem 1.1 (about 1 second; prints PASS, exit code 0)
python check_eta.py --delta "2^-32176926823186696" --x "2^-21451284548791130" --beta0 star --theta "1-2^-40" --mu "2^-30" --eps "12/100000" --gamma "408/100000"

# Theorem 1.2 (prints PASS; (FIN-B) holds)
python check_eta.py --delta "2^-318600000000000000400" --x "2^-212400000000000000200" --beta0 pi8 --theta "1-102/106200000000000000101" --mu "1/100000" --eps "1/110000" --gamma "1/714"

# The second, independently written checker, at its point for Theorem 1.2 (prints PASS)
python reviewer_check.py "eps=2^(-318600000000000000400)" "a=2^(-212400000000000000200)" "Lam=2^(106200000000000000000)" "phi=33/10000000000000000000000" "gam=1/100000" "theta=1/110000" "mu=1/714"
```

On Windows keep the double quotes (`^` is a special character in `cmd`). The `expected_output_*.txt` files contain the outputs we obtained. `check_eta.py --json FILE` writes all balls to a JSON file; `output_theorem_1_1.json` and `output_theorem_1_2.json` are these files at the two points. The option `--prec BITS` sets the working precision (default 256); the verdict at the point of Theorem 1.1 is the same at 64, 128, 256, 512 and 2000 bits. With $N_0-1=32176926823186695$ in place of $N_0$ (same other parameters) the checker prints FAIL.

Optional cross-checks (not part of the proof; need `pip install mpmath sympy`):

```bash
python compare_mpmath.py      # independent float prototype vs check_eta.py; prints COMPARISON OK
python compare_phase5.py      # our chain vs the second checker at the same point (Section 6 of the note)
python algebra/check_phase2_algebra.py    # exact sympy checks of the algebraic identities (Appendix A)
python algebra/check_phase2b_algebra.py
python algebra/check_phase2c_algebra.py
python algebra/check_phase2d_algebra.py
python algebra/check_b1_formulas.py       # exact test of formulas (3), (4), (14) on random systems
```

The sympy scripts check only the algebraic identities between the lines of the proofs. The inequalities are proved in the note.

**Test on real systems** (not part of the proof). Every intermediate inequality of Appendix A holds for every system with $\max_j\lvert S_j\rvert\le q$ and any $q\in(1/2,1/\sqrt2)$, not only for $q$ close to $1/2$. Near-optimal systems with $\max_j\lvert S_j\rvert$ between $0.688$ and $0.707$ exist for $n\ge40$. `stress_test.py` evaluates both sides of 50 of these inequalities, with our explicit constants, on eleven such systems ($40\le n\le1000$) from [turan-power-sum-certificates](https://github.com/Bakeas/turan-power-sum-certificates):

```bash
git clone https://github.com/Bakeas/turan-power-sum-certificates
C=turan-power-sum-certificates
python stress_test.py $C/small_n_tuples.json $C/results/certs/n_300.json $C/results/certs/n_500.json $C/results/certs/n_1000.json
```

It takes about 40 seconds and reports "0 violated" (`expected_output_stress_test.txt`). It cannot see errors that appear only when $q$ is very close to $1/2$.

## Files

| File | Contents |
|---|---|
| `Bakeas_explicit_lower_bound_Turan.pdf`, `.tex` | the note |
| `check_eta.py` | rigorous checker of Theorem 3.1 (python-flint / Arb) |
| `output_theorem_1_1.json`, `output_theorem_1_2.json` | all values computed by `check_eta.py` at the two points (balls) |
| `expected_output_*.txt` | outputs of the checkers |
| `reviewer_check.py` | second rigorous checker, written by a separate AI instance that worked only from Biró's papers, without our derivations or our checker |
| `reviewer_constants.md` | the notes of that instance: its own explicit constants and conditions C1–C11 |
| `compare_mpmath.py`, `compare_phase5.py` | floating-point cross-checks |
| `stress_test.py` | test of the intermediate inequalities on real near-optimal systems |
| `algebra/` | exact computer-algebra checks (sympy) |

**Names used in the code.** The comments in `check_eta.py` refer to "`explicit_constants.md` section 2.27", our working note. Its content is Section 3 (formulas and conditions) and Appendix A (proofs) of the note, with the same numbering of the constants. The second checker uses its own notation: its $\varepsilon$, $a$, $P_t$ are our $\delta$, $x$, $B_t$, and its $C_t$, $D_t$ are our $D_t$, $C_t$; see Section 6 of the note.

**Remark on `reviewer_constants.md`.** Its erratum [E1] (Lemma 1(iii) of Biró's paper printed with ">") came from a low-resolution scan. On the printed page the inequality reads "≥", so [E1] is not an erratum. Its [E2] (in (18), $h_{t-1}$ should be $h_{t+1}$) is a real misprint, found independently by both derivations.

## Related

- V. Bakeas, *Turán's pure power sum problem: certified bounds for finite n and limsup R_n ≤ 0.68637* (upper bounds): [turan-c42-bound](https://github.com/Bakeas/turan-c42-bound).

## AI use

The explicit constants, the proofs, the code and the text were produced with substantial help from an AI model (Claude Opus 5.5, Anthropic), directed by the author, who also ran the checkers. Separate instances of the same model re-derived the constants from Biró's papers without access to our derivations, wrote the second checker, and checked the two main sharpenings (S1 and S8); they found no error. All of these checks were done by instances of one AI model, so they may share a blind spot. No professional mathematician has checked the argument yet.

## License

- Code: MIT (`LICENSE`).
- Text and data: CC BY 4.0.
