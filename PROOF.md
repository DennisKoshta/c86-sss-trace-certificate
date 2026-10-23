# Why the certificate proves C86 ≤ T/m

Notation is as in `FORMAT.md`: μ is a step density with mass m, first moment T and support S. Q_1, …, Q_k
have weights b_j, and f(z) = U_μ(z) − Σ_j b_j log|Q_j(z)| with U_μ(z) = ∫ log|z − y| dμ(y). Put ν = μ/m
(a probability measure) and a_j = b_j/m. Polynomials with b_j = 0 play no role and are dropped.

**Theorem.** If conditions E1–E6 of `FORMAT.md` hold, then C86 ≤ T/m.

The argument is that of Talebizadeh Sardari–Orloski [SO24, Cor. 1.6 and primal problem (11)], resting on
Smith's converse to Fekete's theorem [Smi24], and of the [Pin26] certificate. We give it in full because
two points matter for this certificate: the degree condition must be strict, and the auxiliary
polynomials may have any trace.

## Step 1. f ≥ 0 on all of ℂ

- **Superharmonic off S.** Off S, U_μ is harmonic. Each −b_j log|Q_j| is superharmonic, and +∞ at the
  roots of Q_j. So f is superharmonic on the open set ℂ ∖ S.
- **Lower semicontinuous.** μ has bounded density, so U_μ is continuous on ℂ, and −b_j log|Q_j| is lower
  semicontinuous. Hence liminf_{z→ζ} f(z) ≥ f(ζ) for every ζ ∈ S.
- **Behaviour at ∞.** As z → ∞,

      f(z) = (m − Σ_j b_j deg Q_j) log|z| − Σ_j b_j log|lc Q_j| + o(1),

  which tends to +∞ by the strict inequality E3.

By the minimum principle for superharmonic functions on the Riemann sphere minus S, E5 (f ≥ 0 on S)
gives f ≥ 0 on ℂ. Equivalently,

    (5.17)   U_ν(z) ≥ Σ_j a_j log|Q_j(z)|   for all z ∈ ℂ.

If only Σ b_j deg Q_j = m held, the limit at ∞ would be −Σ_j b_j log|lc Q_j|. That is negative as soon as a
weighted Q_j has |lc| ≥ 2, and then (5.17) fails near ∞. So the strict inequality (or monic Q_j) is needed.
Here every Q_j is monic and Σ b_j deg Q_j / m = 0.83542…, so both readings hold.

## Step 2. ν is arithmetic

ν is arithmetic if ∫ log|P| dν ≥ 0 for every nonzero P ∈ ℤ[x]. This is Smith's Proposition 5.15
[Smi24, arXiv:2111.12660v2, p. 42]. Its hypotheses are:
- the Q_j are distinct, irreducible and primitive (E2);
- a_j > 0 with Σ a_j deg Q_j ≤ 1 (E3);
- (5.17) holds on ℂ (step 1);
- ∫ log|Q_j| dν ≥ 0 (E4).

Smith writes U^μ = −∫ log|z − w| dμ, so his (5.17) U^ν ≤ −Σ a_i log|Q_i| is the inequality above.

**Proof.** ∫ log|·| dν is additive over products, and constants are integers of absolute value ≥ 1, so it
suffices to treat P irreducible of degree n ≥ 1. If P = ±Q_j, use E4. Otherwise let c be the leading
coefficient of P and β_1, …, β_n its roots. Then

    ∫ log|P| dν = log|c| + Σ_i U_ν(β_i)
                ≥ log|c| + Σ_j a_j Σ_i log|Q_j(β_i)|
                = log|c| (1 − Σ_j a_j deg Q_j) + Σ_j a_j log|Res(P, Q_j)|
                ≥ 0.

The middle equality uses Σ_i log|Q_j(β_i)| = log|Res(P, Q_j)| − deg Q_j · log|c|. The last step holds
because Res(P, Q_j) is a nonzero integer (P and Q_j are coprime), |c| ≥ 1 and Σ a_j deg Q_j ≤ 1. ∎

**No condition on the traces of the Q_j appears.** The certificate contains x − 2, x² − 4x + 1,
x² − 4x + 2, x³ − 6x² + 9x − 1 and x³ − 6x² + 9x − 3, whose trace-to-degree ratio is 2, above the bound.
For such P, the inequality ∫ log|P| dν ≥ 0 is one of Smyth's necessary conditions. Putting P into A only
adds that constraint and a weight b_P. [SO24, Cor. 1.6] likewise allows any finite A ⊂ ℤ[x]. Trace
conditions appear in [SO24] only in the uniqueness clause of Thm 1.5, and in [Smi24, Prop. 5.17], about
the optimum. Neither is a hypothesis of the bound.

## Step 3. Smith's Theorem 1.5

[Smi24, arXiv:2111.12660v2, Thm 1.5, p. 4]: let Σ ⊂ ℝ be compact with countably many components and
capacity > 1. Then every arithmetic probability measure supported in Σ is the weak-* limit of the root
measures of distinct monic irreducible integer polynomials R_k whose roots all lie in Σ.

Take Σ = [min S, max(max S, min S + 5)]. This is an interval of length ≥ 5, so its capacity is ≥ 5/4 > 1.
It contains S and lies in (0, ∞) by E1. Here S = [11/200, 519/100], so Σ = [11/200, 519/100], of length
5.135.
- Every root of every R_k lies in Σ ⊂ (0, ∞), so each root of R_k is a totally positive algebraic integer
  of degree deg R_k.
- The function x is continuous on Σ, so tr/deg = ∫ x dμ_{R_k} → ∫ x dν = T/m.
- The R_k are distinct, so these are infinitely many distinct algebraic integers.

Hence C86 ≤ T/m. With E6, C86 ≤ `bound` = 144679917840593158601/79999999999999881200 = 1.80849897300741….

## References

- [Smi24] A. Smith, *Algebraic integers with conjugates in a prescribed distribution*, Ann. of Math. 200
  (2024), 71–122; arXiv:2111.12660 (v2 cited for page numbers).
- [SO24] N. Talebizadeh Sardari, B. J. Orloski, *A quantitative converse of Fekete's theorem*,
  arXiv:2304.10021 (v2, 2024).
- [Pin26] R. Pinho, *A certificate for C86 ≤ 1.8117251*, github.com/rpinho/c86-trace-certificate, commit
  116a58ff99681226ce6299d4a02eeff1e7d3c5c7.
