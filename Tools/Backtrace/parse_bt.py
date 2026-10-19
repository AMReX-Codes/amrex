#!/usr/bin/env python3

import sys
import re
import subprocess
import shutil

if len(sys.argv) < 3:
  print("Usage: %s [exe-file] [backtrace-file]" % sys.argv[0])
  sys.exit(1)

exe_file = sys.argv[1]
bt_file = sys.argv[2]

with open(bt_file, 'rt') as f:
  lines = f.readlines()

frames = []  # (frame, addr)

for l in lines:

  m = re.search(r"libc\.so", l)
  if m:
    continue

  matched = False

  # gnu compiler
  m = re.match(r"\s*(\d+): .*\(\+(0x[\dabcdef]+)\)", l)
  if m:
    matched = True
    frame = m.group(1)
    addr = m.group(2)

  # intel compiler
  if not matched:
    m = re.match(r"\s*(\d+): .*\[(0x[\dabcdef]+)\]", l)
    if m:
      matched = True
      frame = m.group(1)
      addr = m.group(2)

  if matched:
    frames.append((frame, addr))

if not frames:
  sys.exit(0)

# GNU addr2line misreads clang's line-tables-only DWARF 5 output.
# llvm-addr2line handles both gcc and clang output.
tool = shutil.which("llvm-addr2line") or "addr2line"

# One call for all addresses; -a tags each result with its address.
cmd = [tool, "-aCpfie", exe_file] + [addr for _, addr in frames]
out = subprocess.check_output(cmd).decode("utf-8")

results = {}
cur = None
for l in out.splitlines(True):
  m = re.match(r"(0x[\dabcdef]+): (.*\n)", l)
  if m:
    cur = int(m.group(1), 16)
    results[cur] = m.group(2)
  elif cur is not None:
    results[cur] += l

for frame, addr in frames:
  print("%s: %s" % (frame, results.get(int(addr, 16), "")))
