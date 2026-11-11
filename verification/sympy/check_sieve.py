# /// script
# requires-python = ">=3.11"
# dependencies = ["sympy==1.14.0"]
# ///
"""Updates of the multivariate RNA-Sieve (S1, S2)."""

import random

import sympy as sp

from _ledger import N_POINTS, SEED, emit, is_zero, record, spd_from

records = []

# S1. Mean update: stationarity of
# (y - M w)' U^-1 (y - M w) + sum_j C_j (mu_j - m_j)' Sigma_j^-1 (mu_j - m_j)
# at m_j = mu_j + (w_j / C_j) Sigma_j (U + Q)^-1 (y - mu w),
# Q = sum_j (w_j^2 / C_j) Sigma_j.
# (a) scalar case, symbolic.
y, u, w1, w2, c1, c2, s1, s2, mu1, mu2 = sp.symbols(
    "y u w1 w2 C1 C2 s1 s2 mu1 mu2", positive=True
)
m1, m2 = sp.symbols("m1 m2")
objective = (y - m1 * w1 - m2 * w2) ** 2 / u + c1 * (mu1 - m1) ** 2 / s1 + c2 * (
    mu2 - m2
) ** 2 / s2
qq = w1**2 / c1 * s1 + w2**2 / c2 * s2
z = (y - mu1 * w1 - mu2 * w2) / (u + qq)
closed = {m1: mu1 + w1 / c1 * s1 * z, m2: mu2 + w2 / c2 * s2 * z}
ok_scalar = all(is_zero(sp.diff(objective, v).subs(closed)) for v in (m1, m2))

# (b) G = 3, J = 2, exact rational arithmetic.
rng = random.Random(SEED)
G, J = 3, 2
worst = 0
for _ in range(N_POINTS):
    rat = lambda: sp.Rational(rng.randint(1, 9), rng.randint(1, 9))
    U = spd_from([[rat() for _ in range(G)] for _ in range(G)])
    sig = [spd_from([[rat() for _ in range(G)] for _ in range(G)]) for _ in range(J)]
    mu = sp.Matrix(G, J, lambda i, j: rat())
    w = [rat() for _ in range(J)]
    cells = [sp.Integer(rng.randint(5, 500)) for _ in range(J)]
    yv = sp.Matrix([rat() for _ in range(G)])
    Q = sp.zeros(G, G)
    for j in range(J):
        Q += w[j] ** 2 / cells[j] * sig[j]
    zv = (U + Q).inv() * (yv - mu * sp.Matrix(w))
    m_new = sp.Matrix.hstack(
        *[mu[:, j] + w[j] / cells[j] * sig[j] * zv for j in range(J)]
    )
    resid = yv - m_new * sp.Matrix(w)
    for j in range(J):
        grad = -2 * w[j] * U.inv() * resid + 2 * cells[j] * sig[j].inv() * (
            m_new[:, j] - mu[:, j]
        )
        worst = max(worst, max(abs(g) for g in grad))
records.append(
    record(
        "S1",
        "m_j = mu_j + (w_j / C_j) Sigma_j (U + Q)^-1 (y - mu w) is the "
        "stationary point of the mean stage",
        "exact",
        "pass" if ok_scalar and worst == 0 else "fail",
        f"scalar case symbolic; G = 3, J = 2 in exact rationals at {N_POINTS} "
        f"points, largest gradient entry {worst}",
    )
)

# S2. Shrinkage level: off-diagonal terms scale by (1 - l2) / (1 - l1).
l1, l2, d1, d2, rho = sp.symbols("lambda1 lambda2 d1 d2 rho", positive=True)


def shrunk(lam):
    corr = (1 - lam) * sp.Matrix([[1, rho], [rho, 1]]) + lam * sp.eye(2)
    root = sp.diag(sp.sqrt(d1), sp.sqrt(d2))
    return root * corr * root


moved = shrunk(l1)
moved[0, 1] = moved[0, 1] * (1 - l2) / (1 - l1)
moved[1, 0] = moved[1, 0] * (1 - l2) / (1 - l1)
records.append(
    record(
        "S2",
        "Multiplying the off-diagonal terms of the lambda1-shrunk covariance "
        "by (1 - lambda2) / (1 - lambda1) gives the lambda2-shrunk covariance",
        "exact",
        "pass" if is_zero(moved - shrunk(l2)) else "fail",
        "symbolic 2 x 2, estimator D^1/2 {(1 - lambda) R + lambda I} D^1/2",
    )
)

emit(records)
