#!/usr/bin/env python3
"""Independent rigorous verifier for upper-bound certificates of the Schur-Siegel-Smyth trace
constant C86 = lambda_SSS = liminf tr(alpha)/deg(alpha) over totally positive algebraic integers.

Written from the primary sources, sharing no code with the search or with R. Pinho's checker:
  - A. Smith, Algebraic integers with conjugates in a prescribed distribution,
    arXiv:2111.12660v2: Theorem 1.5 (p. 4) and Proposition 5.15 (p. 42);
  - B. Orloski, N. Talebizadeh Sardari, A quantitative converse of Fekete's theorem,
    arXiv:2304.10021v2: primal constraints (11), Corollary 1.6;
  - teorth/optimizationproblems constants/86a.md.

Certificate (JSON object; only these keys are read, every other key is ignored):
  edges      list of N+1 rationals e_0 < ... < e_N
  densities  list of N rationals d_i >= 0 (density of mu on [e_i, e_{i+1}])
  A          list of k integer coefficient lists (constant term first), the polynomials Q_j
  b          list of k rationals b_j >= 0
  bound      rational, must equal int x dmu / mu(R) exactly
A rational is a JSON integer or a string "p", "-p" or "p/q" (decimal digits only). Floats are
rejected wherever a number is read.

Notation. m = mu(R) = sum d_i (e_{i+1} - e_i), T = int x dmu, S = supp mu = union of the cells
with d_i > 0, U(z) = int log|z - y| dmu(y), f(z) = U(z) - sum_j b_j log|Q_j(z)|,
G(u) = u log|u| - u (G(0) = 0; G' = log|u|, G'' = 1/u, G odd), c_k = d_k - d_{k-1}
(d_{-1} = d_N = 0), so that for real x
    U(x) = sum_k c_k G(x - e_k)                      (summation by parts of
                                                      int_{e_i}^{e_{i+1}} log|x-y|dy
                                                      = G(x-e_i) - G(x-e_{i+1})),
and for any z in C off the real axis U(z) = -sum_k c_k Re H(e_k - z), H(w) = w Log w - w
(x - z stays in one open half plane as x runs over R, so Log is analytic along the path).

Theorem used (derived here; conditions checked are V1-V6):
  Put nu = mu/m and a_j = b_j/m. Smith's Prop. 5.15 says: if Sigma is a compact subset of R
  with capacity > 1 and countably many components, nu a probability measure on Sigma, Q_j
  distinct irreducible primitive integer polynomials, a_j > 0, sum a_j deg Q_j <= 1,
      (5.17)  int log|z - w| dnu(w) >= sum a_j log|Q_j(z)|   for ALL z in C,
  and int log|Q_j| dnu >= 0, then int log|P| dnu >= 0 for every nonzero P in Z[x] (proof:
  for P irreducible non-associate to every Q_j, with leading coefficient c and roots beta_i,
  int log|P| dnu = log|c| + sum_i U_nu(beta_i) >= log|c| (1 - sum a_j deg Q_j)
  + sum_j a_j log|Res(P, Q_j)| >= 0). No condition on the traces of the Q_j appears.
  Smith's Thm 1.5 then gives distinct monic irreducible R_k with all roots in Sigma and root
  measures -> nu weak-*. Take Sigma = [min S, max(max S, min S + 5)]: an interval, capacity
  >= 5/4 > 1, inside (0, oo) when min S > 0, so every root of every R_k is a totally positive
  algebraic integer, and tr/deg = int x dmu_{R_k} -> int x dnu = T/m (x is continuous on
  Sigma). Infinitely many distinct such alpha, hence C86 <= T/m. (So the width condition
  max S - min S > 4 is NOT needed: Sigma may be enlarged beyond the support.)
  Reduction of (5.17) to S: off S, U is harmonic and -b_j log|Q_j| superharmonic, so f is
  superharmonic on the open set C \\ S; U is continuous on C (mu has bounded density) and
  -b_j log|Q_j| is lower semicontinuous, so liminf f(z) >= f(zeta) at every zeta in S. At
  infinity f(z) = (m - sum b_j deg Q_j) log|z| - sum b_j log|lc Q_j| + o(1). If
  sum b_j deg Q_j < m this tends to +oo; if equality holds the limit is
  -sum b_j log|lc Q_j|, which is >= 0 only when every Q_j with b_j > 0 has |lc| = 1. In
  those two cases the minimum principle on the Riemann sphere minus S turns f >= 0 on S into
  f >= 0 on C. In the remaining case (equality with some |lc| >= 2) f < 0 near infinity, so
  (5.17) FAILS whatever happens on S: "sum b deg <= m" alone is not enough.

  V1  edges strictly increasing; d_i >= 0; b_j >= 0; m > 0; min S > 0.
  V2  every Q_j (b_j > 0 or not) is a primitive integer polynomial of degree >= 1,
      irreducible over Z (flint factorisation), and no two are associate (Q_i != +-Q_j).
      (Duplicates would in fact be harmless - they merge into one weight - but Prop. 5.15 is
      stated for distinct Q_j, so they are rejected.) Real-rootedness is NOT required.
  V3  sum b_j deg Q_j < m, or = m with |lc Q_j| = 1 whenever b_j > 0 (exact rationals).
  V4  int log|Q_j| dmu = m log|lc Q_j| + sum_{roots r} U(r) >= 0 for every j with b_j > 0
      (arb balls; roots from flint's certified complex root isolation; U(r) closed form).
  V5  f(x) >= 0 for every x in S (continuum check, see below).
  V6  T/m == bound exactly.

V5 method. Each support cell C = [lo, hi] (centre xc, half-width rho) gets a model:
  - far terms (edges with |xc - e_k| >= 8 rho, roots with |xc - r| >= 8 rho) are summed into
    one degree-D Taylor polynomial P_C(x - xc) with arb coefficients
        G(x - e): g_0 = G(u0), g_1 = log|u0|, g_n = (-1)^n / (n (n-1) u0^(n-1)) (n >= 2),
        log|x - r|: l_0 = log|z0|, l_n = Re[(-1)^(n-1) / (n z0^n)] (n >= 1),
    and a rigorous tail bound over |x - xc| <= rho: with t = rho/|u0| <= 1/8,
        |c| |u0| t^(D+1) / ((D+1) D (1-t))    resp.    b t^(D+1) / ((D+1) (1-t));
  - near terms are kept individually.
  A piece J = [a, b] of C (midpoint x0, half-width r, so x = x0 + t, |t| <= r) is certified
  when a lower bound of the form  A + B t + C t^2 - E  is positive on |t| <= r, where
  A, B, C, E collect:
  - P_C: P(x0), P'(x0), P''(x0)/2, and E += max_J |P'''| r^3 / 6 (+ the far tail);
  - a near edge term c G(x - e) with |x0 - e| >= 2 r: its exact order-2 Taylor coefficients
    at x0 and E += |c| |u0| t^3 / (6 (1 - t)), t = r/|u0| <= 1/2;
  - a near edge term touching J (e = a or e = b; e is never inside (a, b) because J lies in
    one cell): x - e has constant sign on J, so c G(x - e) is convex (c (x0 - e) > 0) and
    lies above its tangent at x0, or concave and lies above its chord - both exact linear
    minorants (the slope -oo at e is never evaluated);
  - a near root r of Q_j: order-2 Taylor coefficients at x0 when |x0 - r| >= 2 r, else
    -b_j log|x - r| >= -b_j log max(|a - r|, |b - r|) (|x - r| is maximal at an endpoint).
  The minimum of A + B t + C t^2 over |t| <= r is bounded below by
  A - |B| r + min(0, C_lo) r^2, and by A - B^2 / (4 C_lo) when C_lo > 0 (C_lo = lower end
  of the ball C); the larger is used. Pieces that fail are bisected (exact rational
  midpoints). A piece that still fails at depth MAX_DEPTH is decided by evaluating f at its
  midpoint: f < 0 there gives INVALID with a witness; otherwise the scan of the cell goes on
  (looking for a witness) and the cell ends INVALID (uncertified).
  No decision is ever taken on a float: every comparison is an exact rational comparison or
  an arb ball comparison (which is only True when it holds for every point of the ball).

Output: "VALID <T/m rounded up to 12 decimals>" or "INVALID <reason>"; with --beat V
(decimal or p/q) also "BEATS V" when T/m < V exactly, else "DOES_NOT_BEAT V". --margin X
certifies the stronger f >= X on S instead of f > 0 (a slack report). If that fails, the
certificate is re-checked at level 0: a valid certificate then prints "VALID <bound>" and
"SLACK_BELOW X <reason>" (exit 4), an invalid one prints "INVALID <reason>" (exit 1) as
without --margin. Runtime on Pinho's n = 1200 certificate: about 10 s.

Run with a python that has python-flint (0.9.0 tested), e.g. a venv made by
  uv venv .venv-flint && VIRTUAL_ENV=.venv-flint uv pip install python-flint==0.9.0 mpmath
Exit codes: 0 VALID (and BEATS if --beat given, and f >= X if --margin X given); 1 INVALID;
2 VALID but DOES_NOT_BEAT; 3 usage error (bad arguments, unreadable file, python -O);
4 VALID but SLACK_BELOW the --margin level (takes precedence over 2).
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
import time
from fractions import Fraction

try:
    from flint import acb, arb, arb_poly, ctx, fmpq, fmpz_poly
except ImportError:  # pragma: no cover
    print("usage error: python-flint is required (pip install python-flint==0.9.0)", file=sys.stderr)
    sys.exit(3)

PREC = 128          # bits for all arb work
FAR_FACTOR = 8      # far field: distance >= FAR_FACTOR * rho from the cell centre
DEGREE = 16         # Taylor degree of the far-field model
MAX_DEPTH = 48      # bisection depth limit per cell
MAX_PIECES = 20_000_000
MAX_UNCERTIFIED = 64    # uncertified depth-limit leaves scanned per cell before giving up

_RAT = re.compile(r"-?(0|[1-9][0-9]*)(/[1-9][0-9]*)?")


class CertificateError(ValueError):
    """The certificate is malformed or a condition fails."""


# ----------------------------------------------------------------------------------- parsing

def _rational(v, what: str) -> Fraction:
    if isinstance(v, bool):
        raise CertificateError(f"{what}: boolean where a rational was expected")
    if isinstance(v, int):
        return Fraction(v)
    if isinstance(v, str) and _RAT.fullmatch(v):
        return Fraction(v)
    raise CertificateError(f"{what}: not an exact rational (int or 'p/q' string): {v!r:.60}")


def _int(v, what: str) -> int:
    if isinstance(v, bool) or not isinstance(v, int):
        raise CertificateError(f"{what}: not an integer: {v!r:.60}")
    return v


def _reject_constant(s):
    raise CertificateError(f"non-finite number {s} in JSON")


def parse(text: str):
    try:
        obj = json.loads(text, parse_constant=_reject_constant)
    except json.JSONDecodeError as exc:
        raise CertificateError(f"not valid JSON: {exc}") from None
    if not isinstance(obj, dict):
        raise CertificateError("top level is not a JSON object")
    for key in ("edges", "densities", "A", "b", "bound"):
        if key not in obj:
            raise CertificateError(f"missing key {key!r}")
    for key in ("edges", "densities", "A", "b"):
        if not isinstance(obj[key], list):
            raise CertificateError(f"{key!r} is not a list")
    edges = [_rational(v, f"edges[{i}]") for i, v in enumerate(obj["edges"])]
    dens = [_rational(v, f"densities[{i}]") for i, v in enumerate(obj["densities"])]
    if len(dens) < 1 or len(edges) != len(dens) + 1:
        raise CertificateError(f"need len(edges) == len(densities) + 1 >= 2, got {len(edges)}, {len(dens)}")
    polys = []
    for j, q in enumerate(obj["A"]):
        if not isinstance(q, list) or not q:
            raise CertificateError(f"A[{j}] is not a non-empty list")
        coeffs = [_int(v, f"A[{j}][{i}]") for i, v in enumerate(q)]
        if coeffs[-1] == 0:
            raise CertificateError(f"A[{j}]: leading (last) coefficient is 0")
        polys.append(coeffs)
    weights = [_rational(v, f"b[{j}]") for j, v in enumerate(obj["b"])]
    if len(weights) != len(polys):
        raise CertificateError(f"len(b) = {len(weights)} != len(A) = {len(polys)}")
    bound = _rational(obj["bound"], "bound")
    return edges, dens, polys, weights, bound


# ------------------------------------------------------------------------------ exact checks

def exact_checks(edges, dens, polys, weights, bound):
    """V1, V2, V3, V6. Returns (m, T, support cells, active (poly, b) list, info dict)."""
    N = len(dens)
    for i in range(N):
        if not edges[i] < edges[i + 1]:
            raise CertificateError(f"V1: edges not strictly increasing at index {i}")
        if dens[i] < 0:
            raise CertificateError(f"V1: densities[{i}] < 0")
    for j, w in enumerate(weights):
        if w < 0:
            raise CertificateError(f"V1: b[{j}] < 0")
    m = sum((dens[i] * (edges[i + 1] - edges[i]) for i in range(N)), Fraction(0))
    T = sum((dens[i] * (edges[i + 1] ** 2 - edges[i] ** 2) / 2 for i in range(N)), Fraction(0))
    if m <= 0:
        raise CertificateError("V1: total mass is 0")
    cells = [i for i in range(N) if dens[i] > 0]
    minS, maxS = edges[cells[0]], edges[cells[-1] + 1]
    if not minS > 0:
        raise CertificateError(f"V1: min supp mu = {minS} is not > 0")

    seen = {}
    for j, q in enumerate(polys):
        P = fmpz_poly(q)
        if P.degree() < 1:
            raise CertificateError(f"V2: A[{j}] has degree 0")
        content, factors = P.factor()
        if abs(int(content)) != 1:
            raise CertificateError(f"V2: A[{j}] is not primitive (content {content})")
        if len(factors) != 1 or factors[0][1] != 1:
            raise CertificateError(f"V2: A[{j}] is reducible over Z: {factors}")
        key = tuple(q) if q[-1] > 0 else tuple(-v for v in q)
        if key in seen:
            raise CertificateError(f"V2: A[{j}] is associate to A[{seen[key]}] (duplicate polynomial)")
        seen[key] = j

    active = [(j, polys[j], weights[j]) for j in range(len(polys)) if weights[j] > 0]
    sum_bdeg = sum((w * (len(q) - 1) for _, q, w in active), Fraction(0))
    if sum_bdeg > m:
        raise CertificateError(f"V3: sum b_j deg Q_j = {float(sum_bdeg):.12g} > m = {float(m):.12g}")
    if sum_bdeg == m and any(abs(q[-1]) != 1 for _, q, _ in active):
        raise CertificateError("V3: sum b_j deg Q_j == m with a non-monic Q_j: f(z) -> -sum b log|lc| < 0 at infinity")
    if T / m != bound:
        raise CertificateError(f"V6: bound != T/m (T/m = {float(T / m):.15g})")
    info = {"m": m, "T": T, "sum_bdeg": sum_bdeg, "minS": minS, "maxS": maxS, "cells": len(cells)}
    return m, T, cells, active, info


# ------------------------------------------------------------------------------ arb helpers

def A_(x: Fraction) -> arb:
    return arb(fmpq(x.numerator, x.denominator))


def G(u: Fraction) -> arb:
    """G(u) = u log|u| - u for an exact rational u."""
    if u == 0:
        return arb(0)
    ua = A_(u)
    return ua * A_(abs(u)).log() - ua


def logabs(u: Fraction) -> arb:
    return A_(abs(u)).log()


def roots_of(q):
    """Certified enclosures of all complex roots (Q is squarefree as it is irreducible)."""
    out = []
    for z, mult in fmpz_poly(q).complex_roots():
        if mult != 1:
            raise CertificateError("internal: repeated root of an irreducible polynomial")
        out.append(z)
    if len(out) != len(q) - 1:
        raise CertificateError("internal: root count mismatch")
    return out


def is_real_root(z: acb) -> bool:
    # arb_fmpz_poly_complex_roots returns real roots with exactly zero imaginary part.
    return z.imag.is_zero()


def potential_at_root(z: acb, edges, c) -> arb:
    """U(z) for a root z of some Q_j (real or not), as an arb ball."""
    tot = arb(0)
    if is_real_root(z):
        r = z.real
        for k, ck in c:
            u = r - A_(edges[k])
            if u.is_zero():
                continue                      # G(0) = 0: the root is exactly this edge
            if u.contains(0):
                # |G(u)| <= R log(1/R) + R for |u| <= R < 1/e (u log|u| is monotone there).
                R = arb(u.abs_upper())
                if not R < arb(1) / arb(3):
                    raise CertificateError("internal: root ball too wide")
                g = arb(0, (R * (-R.log()) + R).upper())
            else:
                g = u * abs(u).log() - u
            tot += A_(ck) * g
        return tot
    for k, ck in c:
        w = acb(A_(edges[k])) - z
        tot -= A_(ck) * (w * w.log() - w).real
    return tot


# ---------------------------------------------------------------------- V4: int log|Q| dmu

def check_V4(edges, c, m, active):
    vals = []
    for j, q, w in active:
        val = A_(m) * arb(abs(q[-1])).log()
        for z in roots_of(q):
            val += potential_at_root(z, edges, c)
        if not val > 0:
            if val < 0:
                raise CertificateError(f"V4: int log|Q_{j}| dmu = {val} < 0")
            if not val >= 0:
                raise CertificateError(f"V4: cannot certify int log|Q_{j}| dmu >= 0 (ball {val})")
        vals.append((j, val))
    return vals


# ----------------------------------------------------------------------- V5: f >= 0 on S

class CellModel:
    __slots__ = ("lo", "hi", "xc", "rho", "P", "P1", "P2", "P3", "tail", "near_edges", "near_roots")


def build_cell(lo, hi, edges, c, roots):
    """Far-field Taylor model of f on [lo, hi] (see module docstring)."""
    D = DEGREE
    cm = CellModel()
    cm.lo, cm.hi = lo, hi
    xc = (lo + hi) / 2
    rho = (hi - lo) / 2
    cm.xc, cm.rho = xc, rho
    lim = FAR_FACTOR * rho
    coef = [arb(0) for _ in range(D + 1)]
    S = [arb(0) for _ in range(D + 1)]   # S[n] = sum c_k w_k^(n-1), n >= 2
    tail = arb(0)
    near_edges = []
    for k, ck in c:
        u = xc - edges[k]
        if abs(u) < lim:
            near_edges.append((edges[k], ck))
            continue
        ca = A_(ck)
        ua = A_(u)
        la = logabs(u)
        coef[0] += ca * (ua * la - ua)
        coef[1] += ca * la
        w = 1 / ua
        p = ca * w
        for n in range(2, D + 1):
            S[n] += p
            p *= w
        t = A_(rho / abs(u))
        tail += abs(ca) * abs(ua) * t ** (D + 1) / ((D + 1) * D * (1 - t))
    for n in range(2, D + 1):
        coef[n] += S[n] * ((-1) ** n) / (n * (n - 1))
    tail_a = tail
    near_roots = []
    xca = A_(xc)
    lima = A_(lim)
    rhoa = A_(rho)
    for z, bw, ba in roots:
        if is_real_root(z):
            z0 = xca - z.real
            d = abs(z0)
        else:
            z0 = acb(xca) - z
            d = abs(z0)
        if not d > lima:
            near_roots.append((z, ba))
            continue
        coef[0] -= ba * d.log()
        inv = 1 / z0
        p = inv
        for n in range(1, D + 1):
            term = p.real if isinstance(p, acb) else p
            coef[n] -= ba * term * ((-1) ** (n - 1)) / n
            p *= inv
        t = rhoa / arb(d.lower())
        tail_a += ba * t ** (D + 1) / ((D + 1) * (1 - t))
    cm.P = arb_poly(coef)
    cm.P1 = cm.P.derivative()
    cm.P2 = cm.P1.derivative()
    cm.P3 = cm.P2.derivative()
    cm.tail = arb(tail_a.upper())
    cm.near_edges = near_edges
    cm.near_roots = near_roots
    return cm


def piece_lower_bound(cm: CellModel, a: Fraction, b: Fraction, const: arb) -> arb:
    """Rigorous lower bound of f on [a, b]: an arb ball every point of which is <= inf f on
    [a, b] (see the module docstring)."""
    x0 = (a + b) / 2
    r = (b - a) / 2
    s0 = A_(x0 - cm.xc)
    ra = A_(r)
    A0 = cm.P(s0) + const
    B0 = cm.P1(s0)
    C0 = cm.P2(s0) / 2
    E = cm.tail + arb(cm.P3(arb(s0.mid(), (s0.rad() + ra).upper())).abs_upper()) * ra ** 3 / 6
    Erat = Fraction(0)
    for e, ck in cm.near_edges:
        if a < e < b:
            raise CertificateError("internal: an edge lies inside a piece")
        u0 = x0 - e
        ca = A_(ck)
        if abs(u0) >= 2 * r:
            A0 += ca * G(u0)
            B0 += ca * logabs(u0)
            C0 += ca / (2 * A_(u0))
            t = r / abs(u0)
            Erat += abs(ck) * abs(u0) * t ** 3 / (6 * (1 - t))
        elif ck * u0 > 0:   # convex on [a, b]: tangent at x0
            A0 += ca * G(u0)
            B0 += ca * logabs(u0)
        else:               # concave on [a, b]: chord
            ga, gb = ca * G(a - e), ca * G(b - e)
            A0 += (ga + gb) / 2
            B0 += (gb - ga) / (2 * ra)
    E += A_(Erat)
    x0a = A_(x0)
    for z, ba in cm.near_roots:
        if is_real_root(z):
            z0 = x0a - z.real
        else:
            z0 = acb(x0a) - z
        d = abs(z0)
        if d > 2 * ra:
            inv = 1 / z0
            i1 = inv.real if isinstance(inv, acb) else inv
            i2 = inv * inv
            i2 = i2.real if isinstance(i2, acb) else i2
            A0 -= ba * d.log()
            B0 -= ba * i1
            C0 += ba * i2 / 2
            t = ra / arb(d.lower())
            E += ba * t ** 3 / (3 * (1 - t))
        else:
            if is_real_root(z):
                da, db = abs(A_(a) - z.real), abs(A_(b) - z.real)
            else:
                da, db = abs(acb(A_(a)) - z), abs(acb(A_(b)) - z)
            M = arb(max(da.abs_upper(), db.abs_upper()))
            A0 -= ba * M.log()
    lb = A0 - arb(B0.abs_upper()) * ra - E
    Clo = arb(C0.lower())
    if Clo < 0:
        lb += Clo * ra * ra
    elif Clo > 0:
        lb2 = A0 - B0 * B0 / (4 * Clo) - E
        if lb2.lower() > lb.lower():
            lb = lb2
    return lb


def f_point(x: Fraction, edges, c, roots, const) -> arb:
    """f(x) for a real rational x by the full sums (used only to report witnesses)."""
    val = sum((A_(ck) * G(x - edges[k]) for k, ck in c), arb(0)) + const
    xa = A_(x)
    for z, bw, ba in roots:
        d = abs(xa - z.real) if is_real_root(z) else abs(acb(xa) - z)
        if d.contains(0):
            return arb.pos_inf()
        val -= ba * d.log()
    return val


def check_V5(edges, dens, c, cells, active, log=None, level=arb(0)):
    roots = []
    const = arb(0)
    for j, q, w in active:
        ba = A_(w)
        const -= ba * arb(abs(q[-1])).log()
        for z in roots_of(q):
            roots.append((z, w, ba))
    stats = {"pieces": 0, "max_depth": 0, "min_lb": None, "min_lb_at": None}
    budget = MAX_PIECES
    lvl = "0" if level.is_zero() else level.str(6, radius=False)
    for idx, i in enumerate(cells):
        lo, hi = edges[i], edges[i + 1]
        cm = build_cell(lo, hi, edges, c, roots)
        stack = [(lo, hi, 0)]
        cell_min = None
        uncertified = []
        while stack:
            a, b, dep = stack.pop()
            stats["pieces"] += 1
            budget -= 1
            if budget < 0:
                raise CertificateError("V5: piece budget exhausted")
            lb = piece_lower_bound(cm, a, b, const)
            if lb > level:
                lo_val = lb.lower()
                if cell_min is None or lo_val < cell_min[0]:
                    cell_min = (lo_val, (a + b) / 2)
                if dep > stats["max_depth"]:
                    stats["max_depth"] = dep
                continue
            if dep >= MAX_DEPTH:
                x0 = (a + b) / 2
                fx = f_point(x0, edges, c, roots, const)
                if fx < level:
                    raise CertificateError(f"V5: f(x) < {lvl} at x = {float(x0):.17g} "
                                           f"(exact {x0}); f(x) in {fx.str(10, radius=False)}")
                # Not certified and not refuted: keep scanning this cell for a negative witness
                # (an uncertified leaf is often the edge of a region where f < 0).
                uncertified.append((a, b, fx))
                if len(uncertified) >= MAX_UNCERTIFIED:
                    break
                continue
            mid = (a + b) / 2
            stack.append((mid, b, dep + 1))
            stack.append((a, mid, dep + 1))
        if uncertified:
            a, b, fx = uncertified[0]
            raise CertificateError(f"V5: could not certify f >= {lvl} on [{float(a):.17g}, {float(b):.17g}] "
                                   f"(f(mid) = {fx.str(5, radius=False)}; {len(uncertified)} such pieces "
                                   f"in cell [{float(lo)}, {float(hi)}], no value < {lvl} found)")
        if cell_min is not None and (stats["min_lb"] is None or cell_min[0] < stats["min_lb"]):
            stats["min_lb"], stats["min_lb_at"] = cell_min
        if log and idx % 100 == 0:
            log(f"  V5: cell {idx}/{len(cells)} pieces so far {stats['pieces']}")
    return stats


# ------------------------------------------------------------------------------------- main

def verify(text: str, log=None, level=None):
    ctx.prec = PREC
    edges, dens, polys, weights, bound = parse(text)
    m, T, cells, active, info = exact_checks(edges, dens, polys, weights, bound)
    N = len(dens)
    c = []
    for k in range(N + 1):
        ck = (dens[k] if k < N else 0) - (dens[k - 1] if k > 0 else 0)
        if ck != 0:
            c.append((k, ck))
    if log:
        log(f"  exact checks ok: m = {float(m):.15g}, sum b deg = {float(info['sum_bdeg']):.12g}, "
            f"support [{float(info['minS'])}, {float(info['maxS'])}], {info['cells']} cells, "
            f"{len(active)} weighted Q, {len(c)} density jumps")
    info["V4"] = check_V4(edges, c, m, active)
    if log:
        log("  V4 ok: " + ", ".join(f"Q{j}: {v.str(6, radius=False)}" for j, v in info["V4"]))
    info["V5"] = check_V5(edges, dens, c, cells, active, log, arb(0) if level is None else A_(level))
    return T / m, info


def ceil_dec(x: Fraction, digits: int = 12) -> str:
    s = 10 ** digits
    n = -((-x.numerator * s) // x.denominator)
    sign = "-" if n < 0 else ""
    n = abs(n)
    return f"{sign}{n // s}.{n % s:0{digits}d}"


def main(argv=None) -> int:
    if sys.flags.optimize:
        print("usage error: refusing to run under python -O", file=sys.stderr)
        return 3
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("certificate")
    ap.add_argument("--beat", help="threshold V (decimal or p/q): print BEATS only if bound < V exactly")
    ap.add_argument("-v", "--verbose", action="store_true")
    ap.add_argument("--margin", metavar="X", help="certify the stronger f >= X > 0 on S (slack report; "
                    "VALID then also means min_S f >= X)")
    try:
        args = ap.parse_args(argv)
    except SystemExit:
        return 3
    beat = None
    if args.beat is not None:
        if not re.fullmatch(r"-?[0-9]+(\.[0-9]+)?|-?[0-9]+/[1-9][0-9]*", args.beat):
            print(f"usage error: --beat {args.beat!r} is not a decimal or p/q", file=sys.stderr)
            return 3
        beat = Fraction(args.beat)
    try:
        with open(args.certificate, encoding="utf-8") as fh:
            text = fh.read()
    except OSError as exc:
        print(f"usage error: {exc}", file=sys.stderr)
        return 3
    margin = None
    if args.margin is not None:
        if not re.fullmatch(r"[0-9]+(\.[0-9]+)?|[0-9]+/[1-9][0-9]*", args.margin):
            print(f"usage error: --margin {args.margin!r} is not a decimal or p/q", file=sys.stderr)
            return 3
        margin = Fraction(args.margin)
    log = (lambda s: print(s, file=sys.stderr, flush=True)) if args.verbose else None
    t0 = time.time()
    slack_failure = None
    try:
        value, info = verify(text, log, margin)
    except CertificateError as exc:
        if margin is None or not str(exc).startswith("V5:"):
            print(f"INVALID {exc}")
            return 1
        slack_failure = exc
        if log:
            log(f"  f >= {args.margin} not certified ({exc}); re-checking at level 0")
        try:
            value, info = verify(text, log, None)
        except CertificateError as exc0:
            print(f"INVALID {exc0}")
            return 1
    st = info["V5"]
    if log:
        log(f"  V5 ok: {st['pieces']} pieces, max depth {st['max_depth']}, "
            f"min certified lower bound {st['min_lb'].str(4, radius=False)} near x = {float(st['min_lb_at']):.9g}; "
            f"{time.time() - t0:.1f} s")
    if slack_failure is not None:
        print(f"VALID {ceil_dec(value)}")
        print(f"SLACK_BELOW {args.margin} {slack_failure}")
    else:
        print(f"VALID {ceil_dec(value)}" + (f" (and f >= {args.margin} on S)" if margin else ""))
    code = 4 if slack_failure is not None else 0
    if beat is not None:
        if value < beat:
            print(f"BEATS {args.beat}")
        else:
            print(f"DOES_NOT_BEAT {args.beat}")
            code = code or 2
    return code


if __name__ == "__main__":
    sys.exit(main())
