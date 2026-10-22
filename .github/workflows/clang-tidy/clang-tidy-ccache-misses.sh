#!/usr/bin/env bash
#
# Run clang-tidy on the ccache misses of a CI job, except translation units
# whose hash is in a clang-tidy database that .github/actions/clang-tidy-db/fetch
# put in $RUNNER_TEMP/clang-tidy-db/in.  On success, this run's database is in
# $RUNNER_TEMP/clang-tidy-db/out for .github/actions/clang-tidy-db/save.
#
# usage: clang-tidy-ccache-misses.sh <clang-tidy> <clang-tidy arguments>

set -u
shopt -s nullglob

CT=$1
CT_ARGS=$2
LOG=${CCACHE_LOGFILE:-${GITHUB_WORKSPACE}/ccache.log.txt}
NAME=clang-tidy-db-${GITHUB_WORKFLOW}-${GITHUB_JOB}
DB=${RUNNER_TEMP}/clang-tidy-db
mkdir -p "${DB}/in" "${DB}/out"

"${GITHUB_WORKSPACE}/Tools/C_scripts/mmclt.py" --input "${LOG}" \
    --db-in "${DB}"/in/*/*.txt --db-out "${DB}/out/${NAME}.txt" --db-dir "${DB}/hashes" \
    --clang-tidy="${CT}" --clang-tidy-args="${CT_ARGS}" || exit 1

make -j4 -k -f clang-tidy-ccache-misses.mak
