#!/usr/bin/env python3
"""Run every symbolic check and write the verification ledger.

Usage, from the repository root: python3 verification/run_all.py

Each file `verification/sympy/check_*.py` is run with `uv run --script`
(PEP 723 metadata pins SymPy) and prints its ledger rows as JSON. This
script writes `verification/ledger.json` and
`report/verification-ledger.tex`. Exit status 0 means no row has the
status `fail`. No Lean project is attached to this repository; the
ledger records that.
"""

import datetime
import json
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
CHECKS = sorted((ROOT / "verification" / "sympy").glob("check_*.py"))


def tex(text):
    for old, new in (
        ("\\", "\\textbackslash{}"),
        ("_", "\\_"),
        ("^", "\\^{}"),
        ("&", "\\&"),
        ("%", "\\%"),
        ("#", "\\#"),
    ):
        text = text.replace(old, new)
    return text


def main():
    rows = []
    for check in CHECKS:
        run = subprocess.run(
            ["uv", "run", "--script", str(check)],
            capture_output=True,
            text=True,
            cwd=check.parent,
        )
        if run.returncode != 0:
            rows.append(
                {
                    "id": check.stem,
                    "claim": "check script did not run",
                    "kind": "error",
                    "status": "fail",
                    "detail": run.stderr.strip().splitlines()[-1:],
                    "file": check.name,
                }
            )
            continue
        for row in json.loads(run.stdout.strip().splitlines()[-1]):
            row["file"] = check.name
            rows.append(row)
    ledger = {
        "date": datetime.date.today().isoformat(),
        "engine": "sympy==1.14.0 (mpmath for numeric rows)",
        "lean": "not configured for this repository",
        "rows": rows,
    }
    (ROOT / "verification" / "ledger.json").write_text(
        json.dumps(ledger, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    lines = [
        "% Written by verification/run_all.py; do not edit.",
        "\\begin{longtable}{llp{0.42\\linewidth}lp{0.28\\linewidth}}",
        "\\toprule",
        "Id & Status & Claim & Kind & Detail \\\\",
        "\\midrule",
        "\\endhead",
    ]
    for row in rows:
        lines.append(
            " & ".join(
                tex(str(row[k])) for k in ("id", "status", "claim", "kind", "detail")
            )
            + " \\\\"
        )
    lines += ["\\bottomrule", "\\end{longtable}", ""]
    report = ROOT / "report"
    report.mkdir(exist_ok=True)
    (report / "verification-ledger.tex").write_text("\n".join(lines), encoding="utf-8")
    failed = [row["id"] for row in rows if row["status"] == "fail"]
    for row in rows:
        print(f"{row['id']:>4}  {row['status']:<5} {row['kind']:<18} {row['claim'][:70]}")
    print(f"{len(rows)} rows, {len(failed)} failed.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
