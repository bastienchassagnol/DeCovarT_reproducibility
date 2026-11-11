# /// script
# requires-python = ">=3.11"
# dependencies = ["sympy==1.14.0"]
# ///
"""Covariance of a bulk whose cell-type counts are multinomial (C1-C3).

Notation: J cell types, G genes, p the proportions, M = (m_1 .. m_J) the
mean profiles, Sigma_j the within-type covariances, n cells in the bulk.
"""

import sympy as sp

from _ledger import emit, is_zero, record

J, G = 3, 2
n = sp.Symbol("n", positive=True)
p1, p2 = sp.symbols("p1 p2", positive=True)
p = sp.Matrix([p1, p2, 1 - p1 - p2])
M = sp.Matrix(G, J, lambda g, j: sp.Symbol(f"m{g+1}{j+1}", positive=True))
Sig = [
    sp.Matrix(
        [
            [sp.Symbol(f"s{j+1}_11", positive=True), sp.Symbol(f"s{j+1}_12")],
            [sp.Symbol(f"s{j+1}_12"), sp.Symbol(f"s{j+1}_22", positive=True)],
        ]
    )
    for j in range(J)
]


def total_covariance():
    """S_T(p) as coded in bulk_covariance(total = 1)."""
    out = sp.zeros(G, G)
    for j in range(J):
        out += p[j] * (Sig[j] + M[:, j] * M[:, j].T)
    return out - (M * p) * (M * p).T


records = []

# C1. Derive Cov(N) from the multinomial moment generating function, then
# apply the law of total covariance to y = sum over cells.
t = sp.symbols("t1:4")
mgf = sum(p[j] * sp.exp(t[j]) for j in range(J)) ** n
zero = {s: 0 for s in t}
mean_n = sp.Matrix([sp.diff(mgf, t[j]).subs(zero) for j in range(J)])
second = sp.Matrix(
    J, J, lambda i, j: sp.diff(mgf, t[i], t[j]).subs(zero)
)
cov_n = (second - mean_n * mean_n.T).applyfunc(sp.simplify)
within = sp.zeros(G, G)
for j in range(J):
    within += mean_n[j] * Sig[j]
derived = within + M * cov_n * M.T
ok = is_zero(derived - n * total_covariance())
records.append(
    record(
        "C1",
        "Bulk covariance with multinomial cell-type counts is n S_T(p)",
        "exact",
        "pass" if ok else "fail",
        "law of total covariance with Cov(N) from the multinomial MGF, "
        "G = 2, J = 3; conditional on the within-type moments",
    )
)

# C2. The extra term is the covariance of the type means; rank <= J - 1.
between = sp.zeros(G, G)
centre = M * p
for j in range(J):
    between += p[j] * (M[:, j] - centre) * (M[:, j] - centre).T
plain = sp.zeros(G, G)
for j in range(J):
    plain += p[j] * Sig[j]
ok_form = is_zero(total_covariance() - plain - between)
kernel = (sp.diag(*p) - p * p.T) * sp.ones(J, 1)
ok_rank = is_zero(kernel)
m_num = sp.Matrix(4, 3, [2, 7, 1, 5, 3, 8, 6, 4, 9, 1, 2, 5])
p_num = sp.Matrix([sp.Rational(1, 2), sp.Rational(1, 3), sp.Rational(1, 6)])
rank = (m_num * (sp.diag(*p_num) - p_num * p_num.T) * m_num.T).rank()
records.append(
    record(
        "C2",
        "S_T(p) - sum_j p_j Sigma_j = sum_j p_j (m_j - Mp)(m_j - Mp)'; rank <= J - 1",
        "exact",
        "pass" if ok_form and ok_rank and rank == 2 else "fail",
        f"identity exact; (diag p - pp')1 = 0 exact; rank {rank} for G = 4, J = 3",
    )
)

# C3. Rank-one downdate: determinant and inverse of T - mm'.
a, b, c, u, v = sp.symbols("a b c u v", positive=True)
T = sp.Matrix([[a, b], [b, c]])
m = sp.Matrix([u, v])
q = (m.T * T.inv() * m)[0]
ok_det = is_zero((T - m * m.T).det() - T.det() * (1 - q))
ok_inv = is_zero(
    (T - m * m.T).inv() - (T.inv() + T.inv() * m * m.T * T.inv() / (1 - q))
)
records.append(
    record(
        "C3",
        "det(T - mm') = det(T)(1 - m'T^-1 m) and the matching inverse",
        "exact",
        "pass" if ok_det and ok_inv else "fail",
        "symbolic 2 x 2; S_T is positive definite when m'T^-1 m < 1",
    )
)

emit(records)
