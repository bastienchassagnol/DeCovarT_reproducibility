"""Shared helpers of the SymPy checks: one record per claim, printed as JSON."""

import json
import random

import sympy as sp

SEED = 20261111
N_POINTS = 25


def record(ident, claim, kind, status, detail):
    """One ledger row."""
    return {
        "id": ident,
        "claim": claim,
        "kind": kind,
        "status": status,
        "detail": detail,
    }


def is_zero(expr):
    """Exact test of the rule: simplify, else cancel(together(expand))."""
    expr = sp.sympify(expr)
    if isinstance(expr, sp.MatrixBase):
        return all(is_zero(e) for e in expr)
    if sp.simplify(expr) == 0:
        return True
    return sp.cancel(sp.together(sp.expand(expr))) == 0


def rational_points(symbols, n_points=N_POINTS, low=1, high=9, seed=SEED):
    """Random positive rationals for each symbol, fixed seed."""
    rng = random.Random(seed)
    for _ in range(n_points):
        yield {
            s: sp.Rational(rng.randint(low, high), rng.randint(low, high))
            for s in symbols
        }


def spd_from(rows):
    """A'A + I for a rational matrix A: symmetric positive definite."""
    a = sp.Matrix(rows)
    return a.T * a + sp.eye(a.shape[1])


def emit(records):
    print(json.dumps(records, ensure_ascii=False))
