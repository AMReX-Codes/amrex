#!/usr/bin/env python

import sys

if sys.version_info < (2, 7):
    sys.exit("ERROR: need python 2.7 or later for dep.py")

import argparse
import os
import subprocess


def runcommand(command):
    p = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    out = p.stdout.read()
    return out.strip().decode("ascii")


def get_git_hash(d):
    cwd = os.getcwd()
    os.chdir(d)
    try:
        hash = runcommand("git describe --always --tags --dirty")
    except:
        hash = ""
    os.chdir(cwd)
    return hash


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--git_dirs",
                        help="the directories whose git hashes we should capture",
                        type=str, default="")

    parser.add_argument("--git_names",
                        help="the names corresponding to the directories we hash",
                        type=str, default="")

    # parse and convert to a dictionary
    args = parser.parse_args()

    # git hashes
    if args.git_dirs == "":
        git_dirs = []
    else:
        git_dirs = args.git_dirs.split()

    if args.git_names == "":
        git_names = []
    else:
        git_names = args.git_names.split()

    git_hashes = []
    for d in git_dirs:
        if d and os.path.isdir(d):
            git_hashes.append(get_git_hash(d))
        else:
            git_hashes.append("")

    print("\nsource git hashes:")
    for name, ghash in zip(git_names, git_hashes):
        print("  {:20}: {}".format(name, ghash))
    print("\n")
