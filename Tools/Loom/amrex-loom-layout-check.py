#!/usr/bin/env python3
"""Cross-check launch-tag layouts between the host compiler and the Loom device pass.

The host passes each kernel's launch tag (a struct of lambdas) by copying bytes
at offsets the device pass computed with its own C++ frontend. If the two
frontends lay out a closure differently, the kernel silently reads the wrong
bytes. This script compares, for every kernel in the device-pass manifest, the
leaf fields (offset, size) of the tag against clang's own record layout
(`-Xclang -fdump-record-layouts` from the host compile).

usage: amrex-loom-layout-check.py <clang-layout-dump> <manifest>

stdout: export names that could not be verified against a host layout (one per line)
stderr: a diagnostic per mismatch, and a summary line
exit:   0 always (the driver decides what a mismatch means)
"""
import os
import re
import sys

SCALAR_SIZES = {
    "_Bool": 1, "bool": 1, "char": 1, "signed char": 1, "unsigned char": 1,
    "char8_t": 1, "short": 2, "unsigned short": 2, "char16_t": 2, "int": 4,
    "unsigned int": 4, "char32_t": 4, "wchar_t": 4, "float": 4, "long": 8,
    "unsigned long": 8, "long long": 8, "unsigned long long": 8, "double": 8,
    "long double": 16, "__int128": 16, "unsigned __int128": 16,
}

LINE = re.compile(r"^\s*(\d+)(?::\d+-\d+)?\s\|\s(\s*)(.*)$")
LAMBDA = re.compile(r"\(lambda at ([^)]+)\)")


def norm_loc(loc):
    path, line, col = loc.rsplit(":", 2)
    return f"{os.path.normpath(path)}:{line}:{col}"


def parse_dump(path):
    """Return {closure-location tuple: [layout, ...]} for amrex_loom::Tag records."""
    records = {}
    block = []

    def finish():
        if not block:
            return
        header = block[0]
        m = LINE.match(header)
        if not m or "amrex_loom::Tag<" not in m.group(3):
            return
        if not LAMBDA.search(m.group(3)):
            return
        fields = []  # (offset, depth, text)
        size = None
        for line in block[1:]:
            mm = LINE.match(line)
            if mm:
                fields.append((int(mm.group(1)), len(mm.group(2)) // 2, mm.group(3)))
                continue
            sm = re.search(r"\[sizeof=(\d+)", line)
            if sm:
                size = int(sm.group(1))
        # Key: every closure in the tag, in field pre-order (the device pass
        # emits the same list), e.g. the launch lambda then the user lambda.
        key = tuple(norm_loc(x) for _, _, text in fields for x in LAMBDA.findall(text))
        if key:
            records.setdefault(key, []).append((flatten(fields, size), size))

    with open(path, errors="replace") as f:
        for line in f:
            if line.startswith("*** Dumping AST Record Layout"):
                finish()
                block = []
            elif LINE.match(line) or "[sizeof=" in line:
                block.append(line.rstrip("\n"))
    finish()
    return records


def leaf_size(text):
    """Size of a leaf field from its printed type, or None if not derivable."""
    # Printed as "<type> <name>" (closure captures have no name).
    t = text
    if "*" in t.split("[")[0] or "&" in t:
        return 8, None
    m = re.match(r"^(.*?)\[(\d+)\](?:\s+\S+)?$", t)
    if m:
        elem = m.group(1).strip()
        n = int(m.group(2))
        es = 8 if elem.endswith("*") else SCALAR_SIZES.get(strip_cv(elem))
        return (es, n) if es else (None, n)
    typ = strip_cv(re.sub(r"\s+[A-Za-z_]\w*$", "", t))
    return SCALAR_SIZES.get(typ), None


def strip_cv(t):
    t = re.sub(r"\b(const|volatile)\b", "", t)
    return " ".join(t.split())


def flatten(fields, size):
    """Leaves as (offset, size-or-None); opaque array-of-record spans as ('span', start, end)."""
    leaves = []
    for i, (off, depth, text) in enumerate(fields):
        if i == 0:
            continue  # the record itself
        nxt = fields[i + 1] if i + 1 < len(fields) else None
        if nxt and nxt[1] > depth:
            continue  # container: its children are listed next
        if "(empty)" in text or "(vtable pointer)" in text:
            continue
        es, n = leaf_size(text)
        if n is not None:
            if es is None:
                end = size
                for o2, d2, _ in fields[i + 1:]:
                    if d2 <= depth:
                        end = o2
                        break
                leaves.append(("span", off, end))
            else:
                leaves.extend((off + k * es, es) for k in range(n))
        else:
            leaves.append((off, es))
    return leaves


def parse_manifest(path):
    kernels = []
    cur = None
    with open(path) as f:
        for line in f:
            w = line.split()
            if not w:
                continue
            if w[0] == "kernel":
                cur = {"name": w[1], "size": None, "align": None, "locs": [], "args": []}
                kernels.append(cur)
            elif cur is None:
                continue
            elif w[0] == "size":
                cur["size"], cur["align"] = int(w[1]), int(w[2])
            elif w[0] == "loc":
                cur["locs"].append(norm_loc(w[1]))
            elif w[0] == "arg" and w[1] == "0":
                cur["args"].append((int(w[3]), int(w[4])))
    return kernels


def matches(device, host_leaves, host_size):
    if device["size"] is not None and host_size is not None and device["size"] != host_size:
        return f"sizeof {device['size']} vs host {host_size}"
    args = list(device["args"])
    i = 0
    for leaf in host_leaves:
        if leaf[0] == "span":
            _, start, end = leaf
            while i < len(args) and start <= args[i][0] < end:
                i += 1
            continue
        off, sz = leaf
        if i >= len(args):
            return f"host field at offset {off} has no device argument"
        doff, dsz = args[i]
        if doff != off or (sz is not None and dsz != sz):
            return f"field {i}: device offset {doff} size {dsz}, host offset {off} size {sz}"
        i += 1
    if i != len(args):
        return f"device has {len(args) - i} extra argument(s) from offset {args[i][0]}"
    return None


def split_loc(loc):
    path, line, col = loc.rsplit(":", 2)
    return path, int(line), int(col)


def resolver(records):
    """Map a device-pass lambda position to clang's name for that lambda.

    clang names a lambda written inside macro arguments (AMReX's
    AMREX_LAUNCH_KERNEL) by the macro's expansion point; the device pass
    reports the lambda's spelling. Use clang's exact position when present,
    otherwise the nearest preceding clang lambda position in the same file.
    """
    by_file = {}
    for key in records:
        for loc in key:
            path, line, col = split_loc(loc)
            by_file.setdefault(path, set()).add((line, col))
    def resolve(loc):
        path, line, col = split_loc(loc)
        best = None
        for pos in by_file.get(path, ()):
            if pos <= (line, col) and (best is None or pos > best):
                best = pos
        return f"{path}:{best[0]}:{best[1]}" if best else loc
    return resolve


def main():
    dump, manifest = sys.argv[1], sys.argv[2]
    records = parse_dump(dump)
    resolve = resolver(records)
    checked = unverified = 0
    for k in parse_manifest(manifest):
        key = tuple(resolve(loc) for loc in k["locs"])
        candidates = records.get(key, []) if key else []
        if not candidates:
            unverified += 1
            print(k["name"])
            print(f"layout-check: NO-HOST-INSTANCE (no host layout for {key or 'tag'}): {k['name'][:120]}",
                  file=sys.stderr)
            continue
        # Several instantiations can share one launch site; the tag passes if
        # its layout equals any of them. Report against a same-sized one.
        problems = [matches(k, leaves, size) for leaves, size in candidates]
        checked += 1
        if all(problems):
            same_size = [p for p, (_, size) in zip(problems, candidates) if size == k["size"]]
            print(k["name"])
            if same_size:
                # A host instance of the same size lays its fields out differently.
                print(f"layout-check: MISMATCH {k['name'][:120]}: {same_size[0]}", file=sys.stderr)
            else:
                # No host instance of this size here: the host may not instantiate
                # it in this unit (e.g. extern template), or the size differs.
                print(f"layout-check: NO-HOST-INSTANCE {k['name'][:120]}: sizeof {k['size']}, "
                      f"host sizes {sorted({c[1] for c in candidates})}", file=sys.stderr)
    print(f"layout-check: {checked} checked, {unverified} unverified", file=sys.stderr)


if __name__ == "__main__":
    main()
