#!/usr/bin/env python3
"""Compare a CNS Sod run log with cns-expected.txt.

usage: check-cns.py <run.log> [cns-expected.txt]

Checks the time and dt of every step and the step-20 conserved totals against
the recorded single-source run, within a relative tolerance (default 1e-10),
and reports the totals' distance from the CPU reference. A different GPU model
can reorder floating-point reductions, so agreement to round-off, not bit
identity, is the pass criterion. Exit status 0 means pass.
"""
import os
import re
import sys

RTOL = float(os.environ.get("CNS_RTOL", "1e-10"))


def parse_expected(path):
    steps, totals, cpu = {}, {}, {}
    for line in open(path):
        f = line.split()
        if not f or f[0].startswith("#"):
            continue
        if f[0] == "step":
            steps[int(f[1])] = (float(f[2]), float(f[3]))
        elif f[0] == "cpu_reference":
            cpu[f[1]] = float(f[2])
        else:
            totals[f[0]] = float(f[1])
    return steps, totals, cpu


def parse_run(path):
    steps, totals = {}, {}
    pat = re.compile(r"STEP = (\d+) TIME = (\S+) DT = (\S+)")
    for line in open(path, errors="replace"):
        m = pat.search(line)
        if m:
            steps[int(m.group(1))] = (float(m.group(2)), float(m.group(3)))
        for key, label in (("mass", "Total mass"), ("xmom", "Total x-momentum"),
                           ("energy", "Total energy")):
            if label in line:
                totals[key] = float(line.split()[-1])
    return steps, totals


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    run = sys.argv[1]
    expected = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "cns-expected.txt")
    esteps, etotals, cpu = parse_expected(expected)
    rsteps, rtotals = parse_run(run)
    ok = True
    worst = 0.0
    for n, (t, dt) in sorted(esteps.items()):
        if n not in rsteps:
            print(f"FAIL step {n}: missing from the run log")
            ok = False
            continue
        rt, rdt = rsteps[n]
        worst = max(worst, rel(rt, t), rel(rdt, dt))
    for key, value in etotals.items():
        if key not in rtotals:
            print(f"FAIL {key}: missing from the run log")
            ok = False
            continue
        worst = max(worst, rel(rtotals[key], value))
    identical = all(rsteps.get(n) == v for n, v in esteps.items()) and \
        all(rtotals.get(k) == v for k, v in etotals.items())
    if worst > RTOL:
        ok = False
    print(f"steps checked: {len(esteps)}; max relative difference from the recorded run: {worst:.3e} "
          f"(tolerance {RTOL:g}); bit-identical: {'yes' if identical else 'no'}")
    for key, value in cpu.items():
        if key in rtotals:
            print(f"  {key}: {rtotals[key]!r} (CPU reference {value!r}, relative {rel(rtotals[key], value):.1e})")
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
