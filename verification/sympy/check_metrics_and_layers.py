# /// script
# requires-python = ">=3.11"
# dependencies = ["sympy==1.14.0"]
# ///
"""Aitchison distance, fraction conversion, ILR chart, covariance layers
(M1, M2, F1, Q1)."""

import sympy as sp

from _ledger import emit, is_zero, record

records = []
J = 3

# Helmert basis: column k has 1/sqrt(k(k+1)) on parts 1..k, -k/sqrt(k(k+1))
# on part k + 1.
V = sp.zeros(J, J - 1)
for k in range(1, J):
    for i in range(k):
        V[i, k - 1] = 1 / sp.sqrt(k * (k + 1))
    V[k, k - 1] = -k / sp.sqrt(k * (k + 1))

# M1. Aitchison distance equals the distance of ILR coordinates.
lp = sp.Matrix(sp.symbols("a1:4"))  # log p
lq = sp.Matrix(sp.symbols("b1:4"))  # log q
clr = lambda v: v - sp.ones(J, 1) * sum(v) / J
d_clr = (clr(lp) - clr(lq)).dot(clr(lp) - clr(lq))
d_ilr = (V.T * (lp - lq)).dot(V.T * (lp - lq))
ok = is_zero(d_clr - d_ilr) and is_zero(V.T * V - sp.eye(J - 1))
records.append(
    record(
        "M1",
        "Aitchison distance (clr) equals the Euclidean distance of ILR "
        "coordinates in the Helmert basis",
        "exact",
        "pass" if ok else "fail",
        "symbolic, J = 3; the basis is orthonormal",
    )
)

# M2. Transcript shares to cell proportions.
p1, p2 = sp.symbols("p1 p2", positive=True)
p = [p1, p2, 1 - p1 - p2]
s = sp.symbols("s1:4", positive=True)
q = [p[j] * s[j] / sum(p[k] * s[k] for k in range(J)) for j in range(J)]
back = [q[j] / s[j] / sum(q[k] / s[k] for k in range(J)) for j in range(J)]
records.append(
    record(
        "M2",
        "p_j = (q_j / s_j) / sum_k (q_k / s_k) inverts q_j = p_j s_j / sum_k p_k s_k",
        "exact",
        "pass" if all(is_zero(back[j] - p[j]) for j in range(J)) else "fail",
        "symbolic, J = 3",
    )
)

# F1. Jacobian of the ILR chart p = softmax(V z) and its degeneracy.
z = sp.Matrix(sp.symbols("z1 z2"))
eta = V * z
soft = sp.Matrix([sp.exp(e) for e in eta]) / sum(sp.exp(e) for e in eta)
jac = soft.jacobian(z)
coded = (sp.diag(*soft) - soft * soft.T) * V
ok_jac = is_zero(jac - coded)
ps = sp.Matrix(sp.symbols("p1:4", positive=True))
jp = (sp.diag(*ps) - ps * ps.T) * V
gram = sp.factor(sp.simplify((jp.T * jp).det().subs(ps[2], 1 - ps[0] - ps[1])))
vanishes = all(
    sp.simplify(gram.subs(ps[k], 0)) == 0 for k in range(2)
) and sp.simplify(gram.subs(ps[1], 1 - ps[0])) == 0
records.append(
    record(
        "F1",
        "dp/dz = (diag p - pp')V; det(J'J) vanishes on every face, so the "
        "ILR information J' I_p J is singular there",
        "exact",
        "pass" if ok_jac and vanishes else "fail",
        f"symbolic, J = 3; det(J'J) = {gram}",
    )
)

# Q1. Three covariance layers of the bulk per cell, one gene:
# y = sum_j p_j (mu_j + B_j) + (1/n) sum_j sum_c eps_jc.
n = 12
props = [sp.Rational(1, 4), sp.Rational(3, 4)]
vb = sp.symbols("vB1 vB2", positive=True)
vw = sp.symbols("vW1 vW2", positive=True)
variance = 0
for j, pj in enumerate(props):
    variance += pj**2 * vb[j]  # one replicate effect, coefficient p_j
    cells = int(n * pj)
    variance += cells * sp.Rational(1, n) ** 2 * vw[j]  # independent cells
stated = sum(props[j] ** 2 * vb[j] + props[j] * vw[j] / n for j in range(2))
records.append(
    record(
        "Q1",
        "Cov = sum_j p_j^2 Sigma_B,j + (1/n) sum_j p_j Sigma_W,j for the "
        "replicate-plus-cell hierarchy",
        "exact",
        "pass" if is_zero(variance - stated) else "fail",
        "variance of a linear combination of independent terms, one gene, "
        "n = 12, p = (1/4, 3/4); conditional on independence across cells "
        "and types",
    )
)

emit(records)
