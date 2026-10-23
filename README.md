# A certificate for C<sub>86</sub> ≤ 1.8084990 (Schur–Siegel–Smyth trace constant)

This repository holds an exactly specified measure and auxiliary polynomials, with checkers, proving

    C_86 = λ_SSS ≤ 144679917840593158601 / 79999999999999881200 = 1.80849897300741… < 1.8084990.

Here C<sub>86</sub> is the smallest limit point of tr(α)/deg(α) over totally positive algebraic integers α
(constant 86a of [Tao's optimization-constants table](https://teorth.github.io/optimizationproblems/constants/86a.html)).
The previous best upper bound was 1.8117251, by R. Pinho [Pin26]. Exact margins:
- 1.8117251 − bound = 3.2261269926·10⁻³;
- [Pin26]'s exact value 144938007612277893/79999999999994000 − bound = 3.2261221462·10⁻³.

The gap to the lower bound 1.80203 [OSS25] shrinks from 9.695·10⁻³ to 6.469·10⁻³.

Not refereed.

## Method

The method is that of Talebizadeh Sardari and Orloski [SO24], which rests on Smith's converse to
Fekete's theorem [Smi24], in the certificate format of [Pin26]. Take a step density μ on [11/200, 519/100]
and integer polynomials Q<sub>j</sub> with weights b<sub>j</sub>. Suppose they satisfy:
- f = U<sub>μ</sub> − Σ b<sub>j</sub> log|Q<sub>j</sub>| ≥ 0 on the support;
- Σ b<sub>j</sub> deg Q<sub>j</sub> < μ(ℝ);
- ∫ log|Q<sub>j</sub>| dμ ≥ 0.

Then C<sub>86</sub> ≤ ∫x dμ / μ(ℝ). `PROOF.md` gives the full argument.

**How the polynomials were chosen.** We took every irreducible monic integer polynomial of degree ≤ 7 with
all roots in [0, 5.3] and trace at most twice the degree: 883 polynomials. We solved the discretised
[SO24] primal linear program over all of them at once, with cells of width 1/200, then rounded the
solution to exact rationals with a safety margin. 19 polynomials receive positive weight:
- the eleven of [Pin26];
- x − 2, x² − 4x + 1, x² − 4x + 2 and x³ − 6x² + 9x − 3, of trace-to-degree ratio 2;
- x⁶ − 11x⁵ + 43x⁴ − 73x³ + 53x² − 15x + 1 and x⁶ − 11x⁵ + 43x⁴ − 72x³ + 51x² − 14x + 1;
- x⁷ − 13x⁶ + 64x⁵ − 150x⁴ + 172x³ − 89x² + 18x − 1 and x⁷ − 13x⁶ + 63x⁵ − 143x⁴ + 157x³ − 78x² + 16x − 1.

Polynomials of trace-to-degree ratio 2, above the bound, are allowed: neither Smith's criterion nor
[SO24, Cor. 1.6] places any condition on the traces of the auxiliary polynomials (`PROOF.md`, step 2).

| Q | b<sub>Q</sub> (exact) | ≈ |
|---|---|---|
| x − 1 | 79216914042837/1000000000000000 | 0.079216914043 |
| x | 494739586843/8000000000000 | 0.061842448355 |
| x − 2 | 4041043491479/100000000000000 | 0.040410434915 |
| x² − 3x + 1 | 17462027687139/500000000000000 | 0.034924055374 |
| x² − 4x + 1 | 4580571915059/250000000000000 | 0.018322287660 |
| x² − 4x + 2 | 10103990454883/1000000000000000 | 0.010103990455 |
| x³ − 5x² + 6x − 1 | 17637921744613/500000000000000 | 0.035275843489 |
| x³ − 6x² + 9x − 1 | 110660999533/200000000000000 | 0.000553304998 |
| x³ − 6x² + 9x − 3 | 69448588053/500000000000000 | 0.000138897176 |
| x⁴ − 7x³ + 13x² − 7x + 1 | 17755144029567/1000000000000000 | 0.017755144030 |
| x⁴ − 7x³ + 14x² − 8x + 1 | 17418022936419/1000000000000000 | 0.017418022936 |
| x⁵ − 9x⁴ + 27x³ − 31x² + 12x − 1 | 12331565309643/1000000000000000 | 0.012331565310 |
| x⁵ − 9x⁴ + 28x³ − 35x² + 15x − 1 | 10368009092229/1000000000000000 | 0.010368009092 |
| x⁵ − 9x⁴ + 26x³ − 29x² + 11x − 1 | 4064736833123/500000000000000 | 0.008129473666 |
| x⁵ − 9x⁴ + 27x³ − 32x² + 13x − 1 | 724976585669/125000000000000 | 0.005799812685 |
| x⁶ − 11x⁵ + 43x⁴ − 73x³ + 53x² − 15x + 1 | 3614401969/1000000000000000 | 0.000003614402 |
| x⁶ − 11x⁵ + 43x⁴ − 72x³ + 51x² − 14x + 1 | 73408487/1000000000000000 | 0.000000073408 |
| x⁷ − 13x⁶ + 64x⁵ − 150x⁴ + 172x³ − 89x² + 18x − 1 | 472219232001/62500000000000 | 0.007555507712 |
| x⁷ − 13x⁶ + 63x⁵ − 143x⁴ + 157x³ − 78x² + 16x − 1 | 1521492933501/250000000000000 | 0.006085971734 |

## Contents

| file | what it is |
|---|---|
| `certificate/c86_A19d7_robust_1.8084990.json` | the certificate (format of [Pin26]; keys other than `edges`, `densities`, `A`, `b`, `bound` are metadata) |
| `FORMAT.md` | the format and the conditions E1–E6 a checker must decide |
| `PROOF.md` | why E1–E6 give C<sub>86</sub> ≤ T/m |
| `checkers/c86_trace.py` | independent rigorous verifier (python-flint, mpmath), written from the format and the papers only |
| `checkers/pinho_referee_check_margin.patch` | a one-test change to [Pin26]'s checker that certifies f ≥ 10⁻⁷ instead of f ≥ 0 |
| `checkers/exact_summary.py` | standard-library recomputation of the exact quantities and margins |
| `checkers/requirements.txt` | python-flint 0.9.0, mpmath |
| `expected/` | the outputs of the runs below |
| `SHA256SUMS` | hashes of every file above |

[Pin26]'s checker is not copied here. Fetch it at its pinned commit, as below.

## How to reproduce

In a clone of this repository (Python ≥ 3.10):

```
sha256sum -c SHA256SUMS
python3 -m venv .venv && .venv/bin/pip install -r checkers/requirements.txt
C=certificate/c86_A19d7_robust_1.8084990.json

python3 checkers/exact_summary.py $C                                     # ALL EXACT CHECKS PASS
.venv/bin/python checkers/c86_trace.py $C --beat 1.8117251                # VALID 1.808498973008 / BEATS 1.8117251, exit 0, ~10 s
.venv/bin/python checkers/c86_trace.py $C --beat 1.8117251 --margin 1/10000000   # VALID ... (and f >= 1/10000000 on S), exit 0

# [Pin26]'s checker, unchanged, at its pinned commit
git clone https://github.com/rpinho/c86-trace-certificate pinho
git -C pinho checkout 116a58ff99681226ce6299d4a02eeff1e7d3c5c7
python3 -m venv pinho/.venv && pinho/.venv/bin/pip install -r pinho/requirements.txt
pinho/.venv/bin/python pinho/referee_check.py $C 1                        # ALL_OK True, ~4-5 min on one core
git -C pinho apply ../checkers/pinho_referee_check_margin.patch
pinho/.venv/bin/python pinho/referee_check.py $C 1                        # ALL_OK True: f > 1e-7 on the support
```

The [Pin26] checker writes its report, `c86_A19d7_robust_1.8084990.referee_report.json`, to the current
directory. Our reports and logs are in `expected/`. Its `ALL_OK` includes the clause `bound ≤ 1.8117251`.

What the runs certify:
- E1–E6 hold (both checkers).
- ∫ log|Q| dμ ≥ 1.99·10⁻⁶ for each Q ∈ A.
- f ≥ 10⁻⁷ on the whole support, by both checkers independently. The certificate therefore does not sit on
  the boundary of the conditions.

## Credits

- A. Smith [Smi24], for the converse to Fekete's theorem (Prop. 5.15 and Thm 1.5).
- N. Talebizadeh Sardari and B. J. Orloski [SO24], for the measure construction and the primal problem.
- R. Pinho [Pin26], for the certificate format, the checker and the eleven polynomials, all of which are
  reused here.
- T. Tao and contributors, for the table.

## Disclosure (AI use)

The certificate was found and prepared with Claude (Anthropic, Claude Opus 5.5), directed by Dennis Koshta,
who takes responsibility for it: the linear-programming search over auxiliary polynomial sets, the
certificate, `exact_summary.py`, the patch, and the text. The independent check was done by a separate
Claude agent. That agent wrote `checkers/c86_trace.py` from the format specification and [Smi24]/[SO24]
only, without access to the construction code or to [Pin26]'s checker. The result is not peer-reviewed,
and no human mathematician has checked the argument.

## References

- [Smi24] A. Smith, *Algebraic integers with conjugates in a prescribed distribution*, Ann. of Math. 200
  (2024), 71–122. [arXiv:2111.12660](https://arxiv.org/abs/2111.12660).
- [SO24] N. Talebizadeh Sardari, B. J. Orloski, *A quantitative converse of Fekete's theorem*.
  [arXiv:2304.10021](https://arxiv.org/abs/2304.10021).
- [OSS25] B. J. Orloski, N. Talebizadeh Sardari, A. Smith, *New lower bounds for the Schur–Siegel–Smyth
  trace problem*, Math. Comp. 94 (2025), 2005–2040.
- [Pin26] R. Pinho, *A certificate for C<sub>86</sub> ≤ 1.8117251*,
  [rpinho/c86-trace-certificate](https://github.com/rpinho/c86-trace-certificate), commit
  `116a58ff99681226ce6299d4a02eeff1e7d3c5c7` (MIT).

## License

MIT; see `LICENSE`. [Pin26]'s `referee_check.py`, which the patch modifies, is Copyright (c) 2026 Ricardo
Pinho under the MIT License, and is not included here.
