#!/usr/bin/env python3
"""
Dashboard Generation CLI Script.

Generates the integrated dashboard data payload and HTML artifact for
the RuFaS Brazilian Lactation Equivalence visualization.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Add project root src/ to sys.path
repo_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(repo_root / "src"))

from rufas_lactation_equivalence.dashboard_builder import (  # noqa: E402
    generate_dashboard_payload,
)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build and export the Lactation Equivalence dashboard payload or HTML."
    )
    parser.add_argument(
        "--output-payload",
        type=Path,
        default=None,
        help="Path to write the preprocessed JSON dashboard payload.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress console output.",
    )

    args = parser.parse_args()

    payload = generate_dashboard_payload()

    if not args.quiet:
        summary = payload.get("summary", {})
        ufs_count = len(payload.get("ufs", {}))
        print("=" * 60)
        print("RuFaS Lactation Equivalence Dashboard Payload Builder")
        print("=" * 60)
        print(f"Total UFs processed:     {ufs_count}")
        print(f"RuFaS Macro Regions:     {', '.join(summary.get('regions', []))}")
        print(f"Average Gower Match:     {summary.get('avg_similarity', 0.0):.2%}")
        print(f"Min / Max Similarity:    {summary.get('min_similarity', 0.0):.2%} / {summary.get('max_similarity', 0.0):.2%}")
        print("=" * 60)

    if args.output_payload:
        args.output_payload.parent.mkdir(parents=True, exist_ok=True)
        with open(args.output_payload, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, ensure_ascii=False)
        if not args.quiet:
            print(f"Payload saved to: {args.output_payload}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
