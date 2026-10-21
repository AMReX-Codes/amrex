#!/usr/bin/env python3

"""
A script processing ccache log file to make Make file for Clang-Tidy

This generates a makefile for clang-tidying ccache-cache-missing files. This
could be used to speed up clang-tidy in CIs.

With --db-out, a ccache miss is also skipped if the hash of its translation
unit is in one of the --db-in files.  The hash covers the preprocessed source
with comments, the compile arguments, clang-tidy's version and arguments, and
the --config-file.  When all targets succeed, the makefile writes the hashes
of all translation units in the log to --db-out.

The hits and misses found in the log are checked against `ccache
--print-stats`, so the statistics must be zeroed (ccache -z) before the
build.  A mismatch, e.g. after a change of the log format, is an error.
"""

import argparse, concurrent.futures, hashlib, os, re, shlex, subprocess, sys

def mmclt(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",
                        help="Ccache log file",
                        default="ccache.log.txt")
    parser.add_argument("--identifier",
                        help="Unique identifier for finding compilation line in the log file",
                        default="Src/Base")
    # We assume Src/Base can be used as an identifier to distinguish amrex code
    # from cmake's temporary files like build/CMakeFiles/CMakeScratch/TryCompile-hw3x4m/test_mpi.cpp
    parser.add_argument("--output",
                        help="Make file for clang-tidy",
                        default="clang-tidy-ccache-misses.mak")
    parser.add_argument("--db-in", nargs="*", default=[],
                        help="Files of hashes of translation units that passed clang-tidy")
    parser.add_argument("--db-out",
                        help="File the makefile writes the hashes to when clang-tidy passes")
    parser.add_argument("--db-dir",
                        help="Scratch directory for the hashes (default: DB_OUT.d)")
    parser.add_argument("--clang-tidy", default="clang-tidy",
                        help="Clang-tidy executable (--db-out only)")
    parser.add_argument("--clang-tidy-args", default="",
                        help="Clang-tidy arguments (--db-out only)")
    parser.add_argument("--jobs", type=int, default=os.cpu_count(),
                        help="Parallel preprocessing jobs (--db-out only)")
    parser.add_argument("--ccache", default="ccache",
                        help="Ccache executable for the statistics check; empty skips it (--db-out only)")
    args = parser.parse_args()

    if args.db_out:
        mmclt_db(args)
        return

    with open(args.input, "r") as fin, open(args.output, "w") as fout:

        fout.write("CLANG_TIDY ?= clang-tidy\n")
        fout.write("override CLANG_TIDY_ARGS += --extra-arg=-Wno-unknown-warning-option --extra-arg=-Wno-ignored-optimization-argument --extra-arg-before=--driver-mode=g++\n")
        fout.write("\n")

        fout.write(".SECONDEXPANSION:\n")
        fout.write("clang-tidy: $$(all_targets)\n")
        fout.write("\t@echo SUCCESS\n\n")

        exe_re = re.compile(r" Executing .*? (-.*{}.*) -c(?: \S+)* -o .* (\S*)".format(args.identifier))

        count = 0
        for line in fin.readlines():
            ret_exe_re = exe_re.search(line)
            if ret_exe_re:
                fout.write("target_{}: {}\n".format(count, ret_exe_re.group(2)))
                fout.write("\t$(CLANG_TIDY) $(CLANG_TIDY_ARGS) $< -- {}\n".format
                           (ret_exe_re.group(1)))
                fout.write("\ttouch target_{}\n\n".format(count))
                count = count + 1

        fout.write("all_targets =")
        for i in range(count):
            fout.write(" target_{}".format(i))
        fout.write("\n\n")

        fout.write("clean:\n\t$(RM) $(all_targets)\n\n")


EXTRA_ARGS = ["--extra-arg=-Wno-unknown-warning-option",
              "--extra-arg=-Wno-ignored-optimization-argument",
              "--extra-arg-before=--driver-mode=g++"]
SOURCE_EXTS = (".cpp", ".cxx", ".cc", ".C", ".c")


def read_compilations(logfile, identifier):
    """Compilations in a ccache log: (cwd, compiler, arguments, source, miss),
    and the numbers of hits and misses in the log regardless of identifier."""
    line_re = re.compile(r"^\[\S+ (\d+)\] (.*)$")
    current = {}
    done = []
    with open(logfile, "r", errors="replace") as fin:
        for line in fin:
            m = line_re.match(line.rstrip("\n"))
            if not m:
                continue
            pid, text = m.groups()
            if text.startswith("=== CCACHE") and "STARTED" in text:
                current[pid] = {"cmd": None, "cwd": None, "miss": False}
            elif pid not in current:
                continue
            elif text.startswith("Command line: "):
                current[pid]["cmd"] = text[len("Command line: "):]
            elif text.startswith("Working directory: "):
                current[pid]["cwd"] = text[len("Working directory: "):]
            elif text.startswith("Executing ") and re.search(r" -c( |$)", text):
                current[pid]["miss"] = True
            elif text.startswith("=== CCACHE DONE"):
                done.append(current.pop(pid))
    done.extend(current.values())

    compilations = []
    nhits = nmisses = 0
    for c in done:
        if not c["cmd"] or " -c " not in c["cmd"] + " ":
            continue
        if c["miss"]:
            nmisses += 1
        else:
            nhits += 1
        if identifier not in c["cmd"]:
            continue
        try:
            tokens = shlex.split(c["cmd"])
        except ValueError:
            tokens = c["cmd"].split()
        while tokens and os.path.basename(tokens[0]).startswith("ccache"):
            tokens = tokens[1:]
        parsed = split_command(tokens)
        if parsed:
            compilations.append((c["cwd"], parsed[0], parsed[1], parsed[2], c["miss"]))
    return compilations, nhits, nmisses


def check_ccache_stats(ccache, nhits, nmisses):
    """Abort unless the hits and misses found in the log match ccache's statistics."""
    try:
        out = subprocess.run([ccache, "--print-stats"], stdout=subprocess.PIPE,
                             check=True).stdout.decode(errors="replace")
    except (OSError, subprocess.CalledProcessError) as e:
        sys.exit("mmclt: cannot run {} --print-stats: {}".format(ccache, e))
    stats = {}
    for line in out.splitlines():
        key, _, value = line.partition("\t")
        if value.strip().isdigit():
            stats[key.strip()] = int(value)
    keys = ("direct_cache_hit", "preprocessed_cache_hit", "cache_miss")
    if any(k not in stats for k in keys):
        sys.exit("mmclt: {} --print-stats lacks {}; its format may have changed"
                 .format(ccache, ", ".join(k for k in keys if k not in stats)))
    hits = stats["direct_cache_hit"] + stats["preprocessed_cache_hit"]
    if (nhits, nmisses) != (hits, stats["cache_miss"]):
        sys.exit("mmclt: found {} hits and {} misses in the log, but ccache reports {} and {}."
                 " Was ccache -z run before the build, or has the log format changed?"
                 .format(nhits, nmisses, hits, stats["cache_miss"]))


def split_command(tokens):
    """(compiler, arguments, source), dropping -c and the output and dependency files."""
    if len(tokens) < 2:
        return None
    compiler, rest = tokens[0], tokens[1:]
    src = None
    if "-c" in rest:
        i = rest.index("-c")
        if i + 1 < len(rest) and rest[i+1].endswith(SOURCE_EXTS):
            src = rest[i+1]
    if src is None:
        srcs = [t for t in rest if t.endswith(SOURCE_EXTS)]
        if not srcs:
            return None
        src = srcs[-1]
    out = []
    skip = False
    for t in rest:
        if skip:
            skip = False
        elif t in ("-o", "-MF", "-MT", "-MQ"):
            skip = True
        elif t in ("-c", "-MD", "-MMD", "-MP", "-MG") or t == src:
            pass
        else:
            out.append(t)
    return compiler, out, src


def tu_hash(common, cwd, compiler, cargs, src):
    """Hash of a translation unit, or None if it cannot be preprocessed."""
    try:
        pp = subprocess.run([compiler] + cargs + ["-E", "-C", src], cwd=cwd,
                            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return None
    h = hashlib.sha256(common)
    h.update("\0".join([compiler] + cargs + [src]).encode())
    h.update(b"\0")
    h.update(pp)
    return h.hexdigest()


def mmclt_db(args):
    ct_args = shlex.split(args.clang_tidy_args) + EXTRA_ARGS

    # Everything that affects clang-tidy's result besides the translation unit
    common = hashlib.sha256()
    common.update(b"mmclt-db-1\0")
    try:
        common.update(subprocess.run([args.clang_tidy, "--version"], stdout=subprocess.PIPE,
                                     stderr=subprocess.DEVNULL).stdout)
    except OSError:
        sys.exit("mmclt: cannot run " + args.clang_tidy)
    common.update(b"\0" + "\0".join(ct_args).encode() + b"\0")
    for a in ct_args:
        if a.startswith("--config-file="):
            with open(a[len("--config-file="):], "rb") as f:
                common.update(f.read())
    common = common.digest()

    known = set()
    for fn in args.db_in:
        if os.path.isfile(fn):
            with open(fn) as f:
                known.update(line.strip() for line in f if line.strip())

    compilations, nhits, nmisses = read_compilations(args.input, args.identifier)
    if args.ccache:
        check_ccache_stats(args.ccache, nhits, nmisses)
    if not compilations:
        sys.exit("mmclt: no compilations with {} in {}".format(args.identifier, args.input))
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        hashes = list(ex.map(lambda c: tu_hash(common, *c[:4]), compilations))

    # A hit only counts as passed if the same unit is not also a miss in this
    # log; otherwise the miss decides.
    miss_hashes = {h for (*_, miss), h in zip(compilations, hashes) if miss and h}
    passed = set()      # hashes known to pass: ccache hits and database entries
    targets = []        # (hash, cwd, source, arguments) to clang-tidy
    target_hashes = set()
    nhit = ndb = 0
    for (cwd, compiler, cargs, src, miss), h in zip(compilations, hashes):
        if not miss:
            nhit += 1
            if h and h not in miss_hashes:
                passed.add(h)
        elif h and h in known:
            ndb += 1
            passed.add(h)
        elif h is None or h not in target_hashes:
            target_hashes.add(h)
            targets.append((h, cwd or ".", os.path.join(cwd or ".", src), cargs))

    db_dir = os.path.abspath(args.db_dir or args.db_out + ".d")
    os.makedirs(os.path.dirname(os.path.abspath(args.db_out)), exist_ok=True)
    os.makedirs(db_dir, exist_ok=True)
    for f in os.listdir(db_dir):
        os.remove(os.path.join(db_dir, f))
    for h in passed:
        open(os.path.join(db_dir, h), "w").close()

    def mk(s):  # in a recipe
        return shlex.quote(s).replace("$", "$$")

    def mkname(s):  # a target or prerequisite
        return s.replace("$", "$$").replace(" ", "\\ ")

    names = []
    with open(args.output, "w") as fout:
        fout.write("override CLANG_TIDY := {}\n".format(mk(args.clang_tidy)))
        fout.write("override CLANG_TIDY_ARGS := {}\n\n".format(" ".join(mk(a) for a in ct_args)))
        fout.write(".SECONDEXPANSION:\n")
        fout.write("clang-tidy: $$(all_targets)\n")
        fout.write("\tls {} > {}\n".format(mk(db_dir), mk(args.db_out)))
        fout.write("\t@echo SUCCESS\n\n")
        for i, (h, cwd, src, cargs) in enumerate(targets):
            name = os.path.join(db_dir, h) if h else os.path.abspath("target_{}".format(i))
            names.append(name)
            fout.write("{}: {}\n".format(mkname(name), mkname(src)))
            fout.write("\tcd {} && $(CLANG_TIDY) $(CLANG_TIDY_ARGS) {} -- {}\n".format
                       (mk(cwd), mk(src), " ".join(mk(a) for a in cargs)))
            fout.write("\ttouch $@\n\n")
        fout.write("all_targets = {}\n\n".format(" ".join(mkname(n) for n in names)))
        # Without a hash, always check.
        fout.write(".PHONY: clang-tidy clean {}\n\n".format(
            " ".join(mkname(n) for (h, *_), n in zip(targets, names) if not h)))
        fout.write("clean:\n\t$(RM) $(all_targets)\n\n")

    print("mmclt: {} translation units, {} ccache hits, {} in the database, {} to check"
          .format(len(compilations), nhit, ndb, len(targets)))


if __name__ == "__main__":
    mmclt(sys.argv)
