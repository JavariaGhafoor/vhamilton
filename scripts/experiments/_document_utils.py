"""Shared helpers for document-ordered experiment runners."""
from __future__ import annotations

import csv
import json
import re
import shlex
import subprocess
from pathlib import Path
from typing import Iterable


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_command(out_dir: Path, command: Iterable[str], name: str = "COMMAND.txt") -> None:
    ensure_dir(out_dir)
    (out_dir / name).write_text(shlex.join([str(x) for x in command]) + "\n")


def run_logged(
    command: list[str],
    out_dir: Path,
    *,
    log_name: str = "run.log",
    command_name: str = "COMMAND.txt",
    check: bool = True,
) -> str:
    """Run a command from project root and store its complete stdout/stderr."""
    ensure_dir(out_dir)
    write_command(out_dir, command, command_name)
    proc = subprocess.run(
        command,
        cwd=project_root(),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    text = proc.stdout
    (out_dir / log_name).write_text(text)
    print(text, end="")
    if check and proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, command, output=text)
    return text


def parse_chem_vqe(text: str) -> dict:
    def one(pattern: str, cast=float, default=None):
        m = re.search(pattern, text, flags=re.MULTILINE)
        return default if not m else cast(m.group(1))

    q_match = re.search(r"^\s*best q\s*=\s*(\[[^\n]+\])", text, flags=re.MULTILINE)
    q = json.loads(q_match.group(1)) if q_match else None
    return {
        "best_q": q,
        "vqe_scaled": one(r"^\s*full VQE energy \(scaled\)\s*=\s*([-+0-9.eE]+)"),
        "ground_hk_scaled": one(r"^\s*exact ground of dyadic H_k \(scaled\)\s*=\s*([-+0-9.eE]+)"),
        "variational_gap_scaled": one(r"^\s*variational gap \(scaled\)\s*=\s*([-+0-9.eE]+)"),
        "queries": one(r"^\s*energy evaluations used\s*=\s*([0-9]+)", int),
        "reference_grid": one(r"^\s*exhaustive target grid\s*=\s*([0-9]+)", int),
        "reduction": one(r"^\s*evaluation reduction\s*=\s*([-+0-9.eE]+)x", float),
        "vqe_mh": one(r"^\s*VQE estimate\s*=\s*([-+0-9.eE]+)"),
        "ground_hk_mh": one(r"^\s*exact ground of dyadic H_k\s*=\s*([-+0-9.eE]+)"),
        "physical_ground_mh": one(r"^\s*exact ground of physical H\s*=\s*([-+0-9.eE]+)"),
        "variational_gap_mh": one(r"^\s*variational gap\s*=\s*([-+0-9.eE]+)"),
        "hamilton_shift_mh": one(r"^\s*H discretization shift\s*=\s*([-+0-9.eE]+)"),
        "total_signed_error_mh": one(r"^\s*total signed VQE error\s*=\s*([-+0-9.eE]+)"),
    }


def write_rows_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    ensure_dir(path.parent)
    if not rows:
        path.write_text("")
        return
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            row = dict(row)
            for key, value in row.items():
                if isinstance(value, (list, dict)):
                    row[key] = json.dumps(value, separators=(",", ":"))
            w.writerow(row)
