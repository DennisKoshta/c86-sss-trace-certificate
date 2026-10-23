#!/usr/bin/env python3
"""Exact-arithmetic summary of a C86 upper-bound certificate (standard library only).

Recomputes from the certificate, in rationals: the mass m, the first moment T, the bound T/m
(and checks it equals the certificate's `bound`), the support S, sum_j b_j deg Q_j / m (strict
degree condition), and a compact interval Sigma for Smith's Theorem 1.5. It then prints the
margins against the listed value 1.8117251, against the exact [Pin26] bound, and the remaining
gap to the lower bound 1.80203. It does NOT check the potential inequality or the energies
(conditions E4, E5 in FORMAT.md); use c86_trace.py or the [Pin26] checker for those.

usage: python3 exact_summary.py CERTIFICATE.json        (exit 0 if every check here passes)
"""

import json
import sys
from fractions import Fraction as F

LISTED = F("1.8117251")
PINHO = F(144938007612277893, 79999999999994000)    # bound_exact of [Pin26] cert_n1200_A11.json
LOWER = F("1.80203")


def rat(v):
    if not isinstance(v, (int, str)) or isinstance(v, bool):
        raise ValueError(f"not an exact rational: {v!r}")
    return F(v)


def dec(x, digits=13):
    """x rounded toward zero to `digits` decimals, as a string (exact)."""
    s = "-" if x < 0 else ""
    x = abs(x)
    q = x.numerator * 10**digits // x.denominator
    return f"{s}{q // 10**digits}.{q % 10**digits:0{digits}d}"


def main(path):
    c = json.load(open(path))
    e = [rat(v) for v in c["edges"]]
    d = [rat(v) for v in c["densities"]]
    A = [[int(a) for a in q] for q in c["A"]]
    b = [rat(v) for v in c["b"]]
    ok = True

    def check(name, cond):
        nonlocal ok
        ok &= bool(cond)
        print(f"{'ok  ' if cond else 'FAIL'} {name}")

    check("edges strictly increasing", all(x < y for x, y in zip(e, e[1:])))
    check("len(densities) == len(edges) - 1", len(d) == len(e) - 1)
    check("densities >= 0, weights b >= 0", min(d) >= 0 and min(b) >= 0)
    check("len(b) == len(A)", len(b) == len(A))
    m = sum(di * (e[i + 1] - e[i]) for i, di in enumerate(d))
    T = sum(di * (e[i + 1] ** 2 - e[i] ** 2) / 2 for i, di in enumerate(d))
    cells = [i for i, di in enumerate(d) if di > 0]
    lo, hi = e[cells[0]], e[cells[-1] + 1]
    bound = T / m
    sbd = sum(bj * (len(q) - 1) for q, bj in zip(A, b))
    check("min S > 0", lo > 0)
    check("bound == T/m exactly", bound == rat(c["bound"]))
    check("sum b deg < m (strict)", sbd < m)
    sig_hi = max(hi, lo + 5)
    check("Sigma = [min S, max(max S, min S + 5)] contains S, length > 4, inside (0, oo)",
          lo > 0 and sig_hi >= hi and sig_hi - lo > 4)

    print(f"polynomials {len(A)} (max degree {max(len(q) - 1 for q in A)}), cells {len(d)}, nonzero {len(cells)}")
    print(f"support [{lo}, {hi}]   Sigma [{lo}, {sig_hi}]")
    print(f"m = {m}")
    print(f"T = {T}")
    print(f"bound = T/m = {bound} = {dec(bound)}...")
    print(f"sum b deg / m = {dec(sbd / m, 7)}...")
    print(f"1.8117251 - bound        = {dec(LISTED - bound)}...")
    print(f"[Pin26] exact - bound    = {dec(PINHO - bound)}...")
    print(f"gap to 1.80203: {dec(LISTED - LOWER, 6)} -> {dec(bound - LOWER, 6)}...")
    print(f"bound < 1.8084990: {bound < F('1.8084990')}")
    print("ALL EXACT CHECKS PASS" if ok else "SOME CHECK FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__.strip().splitlines()[-1], file=sys.stderr)
        sys.exit(3)
    sys.exit(main(sys.argv[1]))
