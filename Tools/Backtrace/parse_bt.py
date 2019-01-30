#!/usr/bin/env python

import sys
import re
import subprocess
import shlex

if len(sys.argv) < 3:
  print("Usage: %s [exe-file] [backtrace-file]" % sys.argv[0])
  sys.exit(1)

exe_file = sys.argv[1]
bt_file = sys.argv[2]

with open(bt_file, 'rt') as f:
  lines = f.readlines()

for l in lines:

  matched = False

  # gnu compiler
  m = re.match("\s*(\d+): .*\(\+(0x[\dabcdef]+)\)", l)
  if m:
    matched = True
    frame = m.group(1)
    addr = m.group(2)

  # intel compiler
  if not matched:
    m = re.match("\s*(\d+): \[(0x[\dabcdef]+)\]", l)
    if m:
      matched = True
      frame = m.group(1)
      addr = m.group(2)

  if matched:
    cmd = "addr2line -Cfie %s %s" % (exe_file, addr)
    info = subprocess.check_output(shlex.split(cmd)).decode("utf-8")
    print("%s: %s" % (frame, info))
