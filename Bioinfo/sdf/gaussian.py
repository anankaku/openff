#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import re

# ---- Your desired Gaussian header (WITHOUT charge/multiplicity) ----
NEW_HEADER = """%nprocshared=8
%mem=16GB
# opt b3lyp/6-31g(d,p)

"""

# Regex for charge/multiplicity line, e.g. "0 1", "-1 2"
CHG_MULT_RE = re.compile(r"^\s*-?\d+\s+\d+\s*$")


def replace_header(com_path: Path) -> None:
    """
    Replace everything before the charge/multiplicity line with NEW_HEADER,
    while preserving the original charge/multiplicity line and geometry.
    """
    text = com_path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    chg_idx = None
    for i, line in enumerate(lines):
        if CHG_MULT_RE.match(line):
            chg_idx = i
            break

    if chg_idx is None:
        print(f"[WARNING] No charge/multiplicity line found: {com_path}")
        return

    chg_mult_line = lines[chg_idx].strip()
    geometry_lines = lines[chg_idx + 1 :]

    # Build new file content
    out = []
    out.append(NEW_HEADER.rstrip())
    out.append("")  # blank line between route/header and title block if you want

    # Title line(s): keep simple and consistent
    out.append("QM dihedral opt")
    out.append("")  # blank line required by Gaussian format

    # Preserve original charge/multiplicity
    out.append(chg_mult_line)

    # Append geometry
    out.extend(geometry_lines)

    # Gaussian expects a trailing newline
    com_path.write_text("\n".join(out) + "\n", encoding="utf-8")
    print(f"[OK] Updated: {com_path}")


def main() -> None:
    root = Path(".").resolve()

    # Recursively find all .com files under root
    com_files = sorted(root.rglob("*.com"))

    if not com_files:
        print(f"No .com files found under: {root}")
        return

    for com in com_files:
        replace_header(com)


if __name__ == "__main__":
    main()
