# /// script
# requires-python = ">=3.11"
# dependencies = ["sympy==1.14.0"]
# ///
"""Behaviour of log det of a mixture of covariances (L1, L2)."""

import random

import mpmath as mp
import sympy as sp

from _ledger import N_POINTS, SEED, emit, record

records = []

# L1. Concavity of log det on positive definite matrices, numerically.
mp.mp.dps = 40
rng = random.Random(SEED)
G, J = 6, 3
smallest = None
for _ in range(N_POINTS):
    mats = []
    for _ in range(J):
        a = mp.matrix(G, G)
        for i in range(G):
            for k in range(G):
                a[i, k] = mp.mpf(rng.uniform(-1, 1))
        mats.append(a.T * a + mp.mpf("0.01") * mp.eye(G))
    w = [mp.mpf(rng.uniform(0.05, 1)) for _ in range(J)]
    total = sum(w)
    w = [x / total for x in w]
    mix = mp.zeros(G, G)
    for j in range(J):
        mix += w[j] * mats[j]
    gap = mp.log(mp.det(mix)) - sum(w[j] * mp.log(mp.det(mats[j])) for j in range(J))
    smallest = gap if smallest is None else min(smallest, gap)
records.append(
    record(
        "L1",
        "log det(sum_j p_j Sigma_j) >= sum_j p_j log det Sigma_j",
        "numeric",
        "pass" if smallest >= 0 else "fail",
        f"mpmath, 40 digits, {N_POINTS} random draws, G = 6, J = 3; smallest "
        f"gap {mp.nstr(smallest, 6)}; a numerical check, not a proof",
    )
)

# L2. Two matrices with complementary near-null directions.
eps, p = sp.symbols("epsilon p", positive=True)
A = sp.diag(1, eps)
B = sp.diag(eps, 1)
gap = sp.log((p * A + (1 - p) * B).det()) - (
    p * sp.log(A.det()) + (1 - p) * sp.log(B.det())
)
limit = sp.limit(gap + sp.log(eps), eps, 0, "+")
ok = sp.simplify(limit - sp.log(p * (1 - p))) == 0
first_order = sp.series(gap + sp.log(eps) - sp.log(p * (1 - p)), eps, 0, 2).removeO()
records.append(
    record(
        "L2",
        "For A = diag(1, eps), B = diag(eps, 1): the gap of L1 is "
        "-log(eps) + log(p(1 - p)) + O(eps)",
        "series",
        "pass" if ok else "fail",
        f"expansion in eps to first order; O(eps) term {sp.simplify(first_order)}",
    )
)

emit(records)
