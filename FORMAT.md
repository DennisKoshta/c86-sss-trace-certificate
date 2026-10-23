# Certificate format and the conditions a checker must decide

C86 = λ_SSS is the smallest limit point of tr(α)/deg(α) over totally positive algebraic integers α
(constant 86a of https://teorth.github.io/optimizationproblems/constants/86a.html). An upper-bound
certificate is a JSON file in the format of R. Pinho's certificate (github.com/rpinho/c86-trace-certificate,
commit `116a58ff99681226ce6299d4a02eeff1e7d3c5c7`, file `cert_n1200_A11.json`), so his `referee_check.py`
runs on it unchanged. `PROOF.md` derives the bound from the conditions below.

## Format

A JSON object. A *rational* is a string `"p"`, `"-p"` or `"p/q"` in decimal digits (or a JSON integer).
Only the five keys below are read; every other key (`bound_float`, `report`, ...) is metadata and must
be ignored.

| key | type | meaning |
|---|---|---|
| `edges` | list of N+1 rationals | e_0 < e_1 < … < e_N, the cell edges |
| `densities` | list of N rationals | d_i ≥ 0, the density of μ on [e_i, e_{i+1}] |
| `A` | list of k integer lists | the auxiliary polynomials Q_j, coefficients **constant term first** |
| `b` | list of k rationals | weights b_j ≥ 0, one per Q_j |
| `bound` | rational | the claimed value ∫x dμ / μ(ℝ) |

The measure is μ = Σ_i d_i · 1_{[e_i, e_{i+1}]}(x) dx. Write

- m = μ(ℝ) = Σ_i d_i (e_{i+1} − e_i) and T = ∫x dμ = Σ_i d_i (e_{i+1}² − e_i²)/2;
- S = supp μ = ∪_{i : d_i > 0} [e_i, e_{i+1}];
- U_μ(z) = ∫ log|z − y| dμ(y) and f(z) = U_μ(z) − Σ_j b_j log|Q_j(z)|.

Closed forms, with G(u) = u·log|u| − u and G(0) = 0:

- U_μ(x) = Σ_i d_i [G(x − e_i) − G(x − e_{i+1})] for real x;
- ∫ log|Q| dμ = m·log|lc(Q)| + Σ_{roots r of Q} U_μ(r).

## Conditions

- **E1.** The edges are strictly increasing, every d_i ≥ 0 and b_j ≥ 0, m > 0, and min S > 0.
- **E2.** Each Q_j is a primitive integer polynomial of degree ≥ 1, irreducible over ℤ, and no two are
  associate (Q_i ≠ ±Q_j). Pinho's checker also requires all roots real; the theorem does not.
- **E3.** Σ_j b_j deg Q_j < m, **strictly**. (Equality would also do if every Q_j with b_j > 0 had
  leading coefficient ±1; with equality and some |lc Q_j| ≥ 2, f < 0 near ∞ and the argument fails.)
- **E4.** ∫ log|Q_j| dμ ≥ 0 for every j with b_j > 0.
- **E5.** f(x) ≥ 0 for **every** x ∈ S, not only at sample points. f = +∞ at the roots of the Q_j. At an
  edge where the density jumps, U_μ has a term c·G(x − e_i) whose slope is −∞ there, so f can have narrow
  cusps that uniform sampling misses. A checker needs rigorous enclosures, such as bisection with a
  Taylor or monotonicity lower bound on each piece.
- **E6.** T/m equals `bound` exactly.

The claim is then C86 ≤ T/m (`PROOF.md`). The spread max S − min S needs no check: Smith's Σ only has to
contain S (`PROOF.md`, step 3).

## The certificate in this repository

`certificate/c86_A19d7_robust_1.8084990.json`:
- `edges`: 1197 values 1/50 + k/200 (k = 0, …, 1196), up to 6;
- `densities`: 1196 values, 633 of them nonzero, support S = [11/200, 519/100];
- `A`: 19 monic irreducible polynomials of degree ≤ 7, all real-rooted;
- m = 199999999999999703/200000000000000000;
- `bound` = T/m = 144679917840593158601/79999999999999881200 = 1.80849897300741…;
- Σ b_j deg Q_j / m = 0.83542…

## How the three checkers decide E1–E6

- `checkers/c86_trace.py` (independent verifier) decides E1–E6 (its V1–V6). E4 uses arb balls on
  certified complex root enclosures. E5 uses bisection with an order-2 Taylor lower bound and exact
  tangent or chord minorants for the edge terms; its docstring states the method in full. No decision is
  taken on a float. `--margin X` certifies f ≥ X on S instead of f ≥ 0.
- Pinho's `referee_check.py` @116a58f decides E1–E6 with arb at 128 bits. E5 bisects each support cell
  with a mean-value bound. It demands E4 > 0 and real roots, and its `ALL_OK` also requires
  `bound ≤ 1.8117251`. `checkers/pinho_referee_check_margin.patch` turns its E5 test into f > 1/10000000.
- `checkers/exact_summary.py` (standard library) recomputes E1, E3 and E6, and the margins, in
  rationals.
