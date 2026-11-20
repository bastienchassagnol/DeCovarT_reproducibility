# /// script
# requires-python = ">=3.11"
# dependencies = ["sympy==1.14.0"]
# ///
"""Post hoc contrasts of the type effect (D1-D5), DeCovarT 3.1 #sec-fs-deg-posthoc."""

import sympy as sp

from _ledger import emit, is_zero, rational_points, record

records = []

# D1. Fisher information of eta_j in an NB2 log-link model with offset,
# cell-means design: derived from the log-likelihood of one cell,
# l = x log(m) - (x + 1/phi) log(1 + phi m) + const, m = s exp(eta).
x, s, eta, phi = sp.symbols("x s eta phi", positive=True)
m = s * sp.exp(eta)
ll = x * sp.log(m) - (x + 1 / phi) * sp.log(1 + phi * m)
# Expected information: -E[d2 l / d eta2], linear in x, so replace x by m.
info = sp.simplify(-sp.diff(ll, eta, 2).subs(x, m))
claim = m / (1 + phi * m)
records.append(
    record(
        "D1",
        "I(eta_j) per cell = m/(1 + phi m) for NB2, log link, offset log s (eq-fs-nb-fisher)",
        "exact",
        "pass" if is_zero(info - claim) else "fail",
        "expected information from the NB2 log-likelihood",
    )
)
# Cross terms vanish: a cell of type j does not involve eta_k (cell-means design).

# D2. The J one-versus-average contrasts sum to zero, so all Delta_j <= 0
# forces all Delta_j = 0 (intersection of the A nulls = omnibus null).
ok = True
for J in (3, 4, 5):
    etas = sp.symbols(f"e1:{J + 1}")
    deltas = [etas[j] - sum(etas[k] for k in range(J) if k != j) / (J - 1) for j in range(J)]
    ok &= is_zero(sum(deltas))
records.append(
    record("D2", "sum_j Delta_j = 0 for J = 3, 4, 5 (eq-fs-posthoc-hypotheses)", "exact",
           "pass" if ok else "fail", "symbolic sums")
)

# D3. B implies A: if eta_j - eta_k > 0 for all k != j then Delta_j > 0,
# since Delta_j is the mean of those differences.
J = 4
etas = sp.symbols("e1:5")
d_j = etas[0] - sum(etas[1:]) / (J - 1)
mean_diff = sum(etas[0] - etas[k] for k in range(1, J)) / (J - 1)
records.append(
    record("D3", "Delta_j = mean over k != j of (eta_j - eta_k), so B_j implies A_j", "exact",
           "pass" if is_zero(d_j - mean_diff) else "fail", "J = 4, symbolic")
)

# D4. Counter-example to A implies B: means 10, 8, 1, 1 give Delta > 0 for
# types 1 and 2, and only type 1 above every other type.
mu = [10, 8, 1, 1]
lg = [sp.log(v) for v in mu]
delta = [lg[j] - sum(lg[k] for k in range(4) if k != j) / 3 for j in range(4)]
pos = [bool(sp.N(d) > 0) for d in delta]
sup = [all(lg[j] > lg[k] for k in range(4) if k != j) for j in range(4)]
records.append(
    record("D4", "means 10, 8, 1, 1: A holds for types 1-2, B for type 1 only", "numeric",
           "pass" if pos == [True, True, False, False] and sup == [True, False, False, False] else "fail",
           f"Delta = {[round(float(sp.N(d)), 3) for d in delta]}")
)

# D5. Two-sided B2 implies A2 with the same sign: Delta_j is linear and odd
# in the differences eta_j - eta_k, so all differences < 0 give Delta_j < 0
# (D3 with eta -> -eta), and the two-sided p of a symmetric statistic is
# 2 min(F(t), 1 - F(t)) = 2 F(-|t|).
neg = d_j.subs({e: -e for e in etas}, simultaneous=True)
t = sp.symbols("t", real=True)
F = (1 + sp.erf(t / sp.sqrt(2))) / 2  # standard normal cdf
two_sided = 2 * F.subs(t, -sp.Abs(t))
one_up = 1 - F
ok_p = all(
    is_zero(sp.simplify(two_sided.subs(t, v) - 2 * sp.Min(one_up.subs(t, v), F.subs(t, v))))
    for v in (sp.Rational(-3, 2), sp.Rational(1, 3), sp.Integer(2))
)
records.append(
    record("D5", "Delta_j(-eta) = -Delta_j(eta), so B2 (all differences one sign) implies A2 with that sign; two-sided p = 2 min(upper, lower) one-sided p",
           "exact", "pass" if is_zero(neg + d_j) and ok_p else "fail", "J = 4 symbolic; z at t = -3/2, 1/3, 2")
)

emit(records)
