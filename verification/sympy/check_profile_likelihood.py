# /// script
# requires-python = ">=3.11"
# dependencies = ["sympy==1.14.0"]
# ///
"""Profile likelihood of DeCovarT and its score (P1-P3).

y ~ N(c M p, phi S(p)). Scales c and phi have closed forms; the score in
p is the one coded in decovart_profile_total() and sieve_loglik().
"""

import sympy as sp

from _ledger import N_POINTS, emit, is_zero, rational_points, record

records = []

# P1. Stationary point in (c, phi) for G = 2 and a generic precision W.
c, phi = sp.symbols("c phi", positive=True)
y = sp.Matrix(sp.symbols("y1 y2", positive=True))
m = sp.Matrix(sp.symbols("mu1 mu2", positive=True))
w11, w22 = sp.symbols("w11 w22", positive=True)
w12 = sp.Symbol("w12")
W = sp.Matrix([[w11, w12], [w12, w22]])
G = 2
log_det_s = sp.Symbol("L")
r = y - c * m
loglik = -sp.Rational(G, 2) * sp.log(phi) - log_det_s / 2 - (r.T * W * r)[0] / (2 * phi)
c_hat = (m.T * W * y)[0] / (m.T * W * m)[0]
r_hat = y - c_hat * m
phi_hat = (r_hat.T * W * r_hat)[0] / G
at = {c: c_hat, phi: phi_hat}
ok_c = is_zero(sp.diff(loglik, c).subs(at))
ok_phi = is_zero(sp.diff(loglik, phi).subs(at))
profile = -sp.Rational(G, 2) * sp.log(phi_hat) - log_det_s / 2 - sp.Rational(G, 2)
ok_value = is_zero(loglik.subs(at) - profile)
records.append(
    record(
        "P1",
        "c = m'Wy / m'Wm and phi = r'Wr / G are the stationary scales; "
        "profile log-likelihood = -(G/2) log phi - (1/2) log det S - G/2",
        "exact",
        "pass" if ok_c and ok_phi and ok_value else "fail",
        "symbolic, G = 2, generic symmetric precision",
    )
)


# P2. Score of the profile log-likelihood in p_k, G = 2, J = 2, with and
# without the between-type term; exact rational arithmetic at 25 points.
def score_check(between, tied):
    ps = sp.symbols("p1 p2", positive=True)
    ys = sp.symbols("y1 y2", positive=True)
    ms = sp.symbols("m11 m12 m21 m22", positive=True)
    ss = sp.symbols("a1 b1 c1 a2 b2 c2", positive=True)
    Mm = sp.Matrix([[ms[0], ms[1]], [ms[2], ms[3]]])
    sig = []
    for j in range(2):
        a, b, cc = ss[3 * j : 3 * j + 3]
        # a, a + b^2 ... : positive definite by construction
        sig.append(sp.Matrix([[a + b**2, b * cc], [b * cc, a + cc**2]]))
    pv = sp.Matrix(ps)
    yv = sp.Matrix(ys)
    total = (ps[0] + ps[1]) if tied else 1

    def cov(pvec):
        out = sp.zeros(2, 2)
        for j in range(2):
            out += pvec[j] * sig[j]
            if between:
                out += pvec[j] * Mm[:, j] * Mm[:, j].T
        if between:
            tot = (pvec[0] + pvec[1]) if tied else 1
            out -= (Mm * pvec) * (Mm * pvec).T / tot
        return out

    S = cov(pv)
    Sinv = S.inv()
    mean = Mm * pv
    if tied:
        # RNA-Sieve scale: mean M w, covariance U(w); no free scale.
        resid = yv - mean
        value = -(resid.T * Sinv * resid)[0] / 2 - sp.log(S.det()) / 2
        a_vec = Sinv * resid
        lead = [(Mm[:, k].T * a_vec)[0] for k in range(2)]
        scale = 1
    else:
        c_hat = (mean.T * Sinv * yv)[0] / (mean.T * Sinv * mean)[0]
        resid = yv - c_hat * mean
        phi_hat = (resid.T * Sinv * resid)[0] / 2
        value = -sp.log(phi_hat) - sp.log(S.det()) / 2 - 1
        a_vec = Sinv * resid
        lead = [c_hat / phi_hat * (Mm[:, k].T * a_vec)[0] for k in range(2)]
        scale = phi_hat
    worst = 0
    for point in rational_points(list(ps) + list(ys) + list(ms) + list(ss)):
        for k in range(2):
            d_s = sig[k]
            if between:
                d_s = d_s + Mm[:, k] * Mm[:, k].T - (
                    Mm[:, k] * mean.T + mean * Mm[:, k].T
                ) / total
                if tied:
                    d_s = d_s + mean * mean.T / total**2
            coded = (
                lead[k]
                + (a_vec.T * d_s * a_vec)[0] / (2 * scale)
                - (Sinv * d_s).trace() / 2
            )
            derived = sp.diff(value, ps[k])
            gap = sp.nsimplify((coded - derived).subs(point))
            if gap != 0:
                worst = max(worst, abs(sp.N(gap)))
    return worst


details = []
status = "pass"
for between in (False, True):
    for tied in (False, True):
        worst = score_check(between, tied)
        details.append(
            f"between={between}, tied={tied}: max |gap| {worst}"
        )
        if worst != 0:
            status = "fail"
records.append(
    record(
        "P2",
        "Score in p_k: (c/phi) m_k'a + a'(dS_k)a / (2 phi) - tr(S^-1 dS_k) / 2, "
        "a = S^-1 r; and its tied-scale form",
        "rational-sampling",
        status,
        f"exact rational arithmetic at {N_POINTS} points, G = 2, J = 2; "
        + "; ".join(details),
    )
)

# P3. With phi fixed at its profiled value q / G, the gradient of the
# fixed-scale quadratic term equals that of the profiled term.
theta = sp.Symbol("theta")
q = sp.Function("q", positive=True)(theta)
Gs = sp.Symbol("G", positive=True)
profiled = sp.diff(-Gs / 2 * sp.log(q), theta)
fixed = sp.diff(-q / (2 * sp.Symbol("phi0", positive=True)), theta).subs(
    sp.Symbol("phi0", positive=True), q / Gs
)
records.append(
    record(
        "P3",
        "d/dtheta[-(G/2) log q] equals d/dtheta[-q / (2 phi)] at phi = q / G: "
        "profiling divides the gradient of the quadratic form by phi-hat",
        "exact",
        "pass" if is_zero(profiled - fixed) else "fail",
        "symbolic, generic q(theta)",
    )
)

emit(records)
