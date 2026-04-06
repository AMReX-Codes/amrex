#!/usr/bin/env python3

from __future__ import annotations

import argparse
import pathlib
import re
import sys
from collections import OrderedDict
from typing import Tuple

HEADING_RE = re.compile(r"^#\s+[0-9]{2}\.[0-9]{2}\s*$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fail if CHANGES.md contains duplicate '# YY.MM' release headings."
        )
    )
    parser.add_argument(
        "changes_file",
        nargs="?",
        default="CHANGES.md",
        help="Path to the changelog (default: %(default)s)",
    )
    return parser.parse_args()


def load_headings(changelog: pathlib.Path) -> OrderedDict[str, list[int]]:
    duplicates: OrderedDict[str, list[int]] = OrderedDict()
    with changelog.open("r", encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, start=1):
            if HEADING_RE.match(line.rstrip("\n")):
                duplicates.setdefault(line.strip(), []).append(lineno)
    return duplicates


def parse_release_value(heading: str) -> Tuple[int, int]:
    # heading is guaranteed to look like "# YY.MM"
    match = re.search(r"(\d{2})\.(\d{2})", heading)
    if not match:
        return (0, 0)
    year = int(match.group(1))
    month = int(match.group(2))
    return year, month


def main() -> int:
    args = parse_args()
    changelog = pathlib.Path(args.changes_file)

    if not changelog.exists():
        print(f"error: '{changelog}' does not exist", file=sys.stderr)
        return 1

    headings = load_headings(changelog)
    dup_list = [
        (heading, lines) for heading, lines in headings.items() if len(lines) > 1
    ]

    if dup_list:
        print(f"Duplicate release headings detected in {changelog}:")
        for heading, lines in dup_list:
            joined = ", ".join(str(num) for num in lines)
            print(f"  {heading} appears {len(lines)} times (lines: {joined})")
        return 1

    releases: list[Tuple[str, Tuple[int, int]]] = []
    month_errors: list[str] = []
    for heading, lines in headings.items():
        value = parse_release_value(heading)
        releases.append((heading, value))
        if not 1 <= value[1] <= 12:
            joined = ", ".join(str(num) for num in lines)
            month_errors.append(f"{heading} (lines: {joined})")

    if month_errors:
        print(
            f"Invalid month values in {changelog}; expected 01-12 after the dot:"
        )
        for entry in month_errors:
            print(f"  {entry}")
        return 1

    previous_value: Tuple[int, int] | None = None
    previous_heading: str | None = None
    for heading, value in releases:
        if previous_value is not None and value >= previous_value:
            print(
                "Release headings must be in strictly descending order "
                "(latest first)."
            )
            print(f"  Found {heading} after {previous_heading}.")
            return 1
        previous_value = value
        previous_heading = heading

    print(f"No duplicate YY.MM headings found in {changelog}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
