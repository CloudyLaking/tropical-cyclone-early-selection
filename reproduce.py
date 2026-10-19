"""Verify frozen inputs and redraw manuscript figures without an archive download."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def verify_inputs():
    """Check every distributed data/geometry file using the standard library."""
    manifest = json.loads((ROOT / "data/manifest.json").read_text(encoding="utf-8"))
    failures = []
    for relative, expected in manifest["files"].items():
        path = ROOT / relative
        if not path.is_file():
            failures.append(relative + " (missing)")
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if path.stat().st_size != expected["bytes"] or digest != expected["sha256"]:
            failures.append(relative + " (changed)")
    if failures:
        raise SystemExit("Input verification failed:\n" + "\n".join(failures))
    print(f"Verified {len(manifest['files'])} input files.", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--main", action="store_true", help="Draw the five main figures (default)."
    )
    group.add_argument(
        "--supplement", action="store_true", help="Draw supporting Figures S1-S15."
    )
    group.add_argument("--all", action="store_true", help="Draw all twenty figures.")
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Only verify inputs unless a drawing option is also given.",
    )
    args = parser.parse_args()
    verify_inputs()
    if args.verify and not (args.main or args.supplement or args.all):
        return
    (ROOT / "results/checks").mkdir(parents=True, exist_ok=True)
    scripts = []
    expected = []
    if args.all or not args.supplement:
        scripts += [f"figure{i}.py" for i in range(1, 6)]
        expected += [f"Figure{i}" for i in range(1, 6)]
    if args.all or args.supplement:
        scripts += ["supplementary.py", "spatial.py"]
        expected += [f"FigureS{i:02d}" for i in range(1, 16)]
    for script in scripts:
        print(f"Rendering {script}", flush=True)
        subprocess.run(
            [sys.executable, str(ROOT / "src" / script)], cwd=ROOT, check=True
        )
    missing = [
        name + ext
        for name in expected
        for ext in (".pdf", ".svg", ".png")
        if not (ROOT / "results" / (name + ext)).is_file()
    ]
    if missing:
        raise SystemExit("Missing outputs: " + ", ".join(missing))
    print(f"Finished {len(expected)} figures. Outputs: {ROOT / 'results'}")


if __name__ == "__main__":
    main()
