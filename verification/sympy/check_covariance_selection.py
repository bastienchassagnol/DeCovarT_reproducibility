# /// script
# requires-python = ">=3.11"
# dependencies = ["sympy==1.14.0"]
# ///
"""Shrinkage path, factor form, restricted likelihood, linked scales,
conditioning along the shrinkage path (R1 to R6)."""

import sympy as sp

from _ledger import emit, is_zero, record

records = []
lam, ell = sp.symbols("lambda ell", positive=True)

# R1. K = (1 - lambda) M + lambda I keeps the eigenvectors of M, with
# eigenvalues (1 - lambda) ell + lambda. Generic symmetric 2 x 2 M.
a, b, c = sp.symbols("a b c", real=True)
m = sp.Matrix([[a, b], [b, c]])
k = (1 - lam) * m + lam * sp.eye(2)
ok = all(
    is_zero(k * vec - ((1 - lam) * val + lam) * vec)
    for val, _, vecs in m.eigenvects()
    for vec in vecs
)
records.append(
    record(
        "R1",
        "(1 - lambda) M + lambda I has the eigenvectors of M and eigenvalues "
        "(1 - lambda) ell + lambda",
        "exact",
        "pass" if ok else "fail",
        "generic symmetric 2 x 2 M, eigenpairs from sympy",
    )
)

# R2. Rank-deficient M = ell v v' in three dimensions: determinant and
# quadratic form of K through the range and the null space.
v = sp.Matrix([1, 2, 2]) / 3
z = sp.Matrix(sp.symbols("z1 z2 z3", real=True))
k = (1 - lam) * ell * v * v.T + lam * sp.eye(3)
e = (1 - lam) * ell + lam
proj = (v.T * z)[0] ** 2
ok = is_zero(k.det() - e * lam**2) and is_zero(
    (z.T * k.inv() * z)[0] - (proj / e + ((z.T * z)[0] - proj) / lam)
)
records.append(
    record(
        "R2",
        "With M of rank r in p dimensions, det K = prod(e_i) lambda^(p - r) "
        "and z' K^-1 z = sum (v_i' z)^2 / e_i + (|z|^2 - sum (v_i' z)^2) / lambda",
        "exact",
        "pass" if ok else "fail",
        "p = 3, r = 1, v = (1, 2, 2) / 3, symbolic ell, lambda and z",
    )
)

# R3. Factor form B B' + Psi: Woodbury inverse and determinant lemma.
bvec = sp.Matrix(sp.symbols("b1 b2 b3", real=True))
psi = sp.diag(*sp.symbols("psi1 psi2 psi3", positive=True))
sigma = bvec * bvec.T + psi
mcap = 1 + (bvec.T * psi.inv() * bvec)[0]
woodbury = psi.inv() - psi.inv() * bvec * bvec.T * psi.inv() / mcap
ok = is_zero(sigma * woodbury - sp.eye(3)) and is_zero(sigma.det() - psi.det() * mcap)
records.append(
    record(
        "R3",
        "(B B' + Psi)^-1 = Psi^-1 - Psi^-1 B (I + B' Psi^-1 B)^-1 B' Psi^-1 "
        "and det(B B' + Psi) = det(Psi) det(I + B' Psi^-1 B)",
        "exact",
        "pass" if ok else "fail",
        "three genes, one factor, symbolic loadings and uniquenesses",
    )
)

# R4. Restricted likelihood. y ~ N(x beta, phi V), beta integrated out
# with a flat prior: n = 2 observations, one coefficient, V diagonal.
y1, y2, beta = sp.symbols("y1 y2 beta", real=True)
x1, x2 = sp.symbols("x1 x2", real=True, nonzero=True)
v1, v2, phi = sp.symbols("v1 v2 phi", positive=True)
density = sp.exp(
    -((y1 - x1 * beta) ** 2 / v1 + (y2 - x2 * beta) ** 2 / v2) / (2 * phi)
) / (2 * sp.pi * phi * sp.sqrt(v1 * v2))
marginal = sp.integrate(density, (beta, -sp.oo, sp.oo))
info = x1**2 / v1 + x2**2 / v2
beta_hat = (x1 * y1 / v1 + x2 * y2 / v2) / info
rss = (y1 - x1 * beta_hat) ** 2 / v1 + (y2 - x2 * beta_hat) ** 2 / v2
stated = -sp.Rational(1, 2) * (
    sp.log(phi) + sp.log(v1 * v2) + sp.log(info) + rss / phi
)
difference = sp.simplify(sp.expand_log(sp.log(marginal), force=True) - stated)
constant = not difference.free_symbols
# Profile in phi: maximiser rss / (n - p), with n - p = 1 here.
stationary = is_zero(sp.diff(stated, phi).subs(phi, rss))
records.append(
    record(
        "R4",
        "Restricted log-likelihood -1/2 [(n - p) log phi + log det V "
        "+ log det(X' V^-1 X) + r' V^-1 r / phi], maximised at "
        "phi = r' V^-1 r / (n - p)",
        "exact",
        "pass" if constant and stationary else "fail",
        "n = 2, p = 1, diagonal V: Gaussian integral over beta by sympy; "
        "conditional on V not depending on beta",
    )
)

# R5. Linked scales: phi = a c (c = a n, phi = a^2 n) with a fixed. The
# log-likelihood in c is -G/2 log(a c) - (Q - 2 c N + c^2 D) / (2 a c);
# N drops out of its derivative, whose positive root is
# c = (-G a + sqrt(G^2 a^2 + 4 D Q)) / (2 D), and the second derivative
# is negative there.
G_, a_, D_, Q_, c_ = sp.symbols("G a D Q c", positive=True)
N_ = sp.symbols("N", real=True)
ell_c = -G_ / 2 * sp.log(a_ * c_) - (Q_ - 2 * c_ * N_ + c_**2 * D_) / (2 * a_ * c_)
score = sp.simplify(sp.diff(ell_c, c_))
root = (-G_ * a_ + sp.sqrt(G_**2 * a_**2 + 4 * D_ * Q_)) / (2 * D_)
no_n = N_ not in score.free_symbols
at_root = is_zero(score.subs(c_, root))
# sqrt(G^2 a^2 + 4 D Q) > G a because 4 D Q > 0: check the squared gap.
gap = sp.expand((G_**2 * a_**2 + 4 * D_ * Q_) - (G_ * a_) ** 2)
concave = sp.simplify(
    sp.diff(ell_c, c_, 2).subs(c_, root) * (2 * a_ * root**3)
)
concave_ok = is_zero(concave + Q_ + D_ * root**2)
records.append(
    record(
        "R5",
        "Linked scales phi = a c, a fixed: the score in c does not contain "
        "N; its root c = (-G a + sqrt(G^2 a^2 + 4 D Q)) / (2 D) is positive "
        "and a maximum",
        "exact",
        "pass" if no_n and at_root and gap == 4 * D_ * Q_ and concave_ok else "fail",
        "symbolic in G, a, D, Q > 0 and real N; positivity because "
        "G^2 a^2 + 4 D Q - (G a)^2 = 4 D Q > 0; second derivative times "
        "2 a c^3 equals G a c - 2 Q = -(Q + D c^2) < 0 at the root",
    )
)

# R6. Conditioning along the shrinkage path of a correlation matrix
# (docs/04, regularisation callout; scripts/04_05_vis_regularisation.R).
# With the eigenvalues of R_lambda = (1 - lambda) R + lambda I given by R1:
# (a) their sum stays G, the trace of a correlation matrix;
# (b) kappa(lambda) = ((1 - lambda) l_max + lambda) / ((1 - lambda) l_min + lambda)
#     has derivative (l_min - l_max) / ((1 - lambda) l_min + lambda)^2 <= 0;
# (c) kappa(1) = 1 and the effective rank (sum l)^2 / sum l^2 is G at lambda = 1.
l1, l2, l3 = sp.symbols("l1 l2 l3", nonnegative=True)
spectrum = [l1, l2, 3 - l1 - l2]  # G = 3, trace 3
shrunk = [(1 - lam) * v + lam for v in spectrum]
trace_ok = is_zero(sum(shrunk) - 3)
lmax, lmin = sp.symbols("l_max l_min", nonnegative=True)
kappa = ((1 - lam) * lmax + lam) / ((1 - lam) * lmin + lam)
dkappa_ok = is_zero(sp.diff(kappa, lam) - (lmin - lmax) / ((1 - lam) * lmin + lam) ** 2)
at_one = is_zero(kappa.subs(lam, 1) - 1) and is_zero(
    (sum(shrunk) ** 2 / sum(v**2 for v in shrunk)).subs(lam, 1) - 3
)
records.append(
    record(
        "R6",
        "Shrinkage path of a correlation matrix: eigenvalue sum stays G; "
        "kappa(lambda) is non-increasing, d kappa / d lambda = "
        "(l_min - l_max) / ((1 - lambda) l_min + lambda)^2; kappa(1) = 1 and "
        "effective rank (sum l)^2 / sum l^2 = G at lambda = 1",
        "exact",
        "pass" if trace_ok and dkappa_ok and at_one else "fail",
        "symbolic in lambda and the eigenvalues; trace and effective rank "
        "checked with G = 3 (eigenvalues l1, l2, 3 - l1 - l2); eigenvalue "
        "map from R1",
    )
)

emit(records)
