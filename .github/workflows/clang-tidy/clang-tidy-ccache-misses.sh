#!/usr/bin/env bash
#
# Run clang-tidy on the ccache misses of a CI job, except translation units
# whose hash is in the clang-tidy database of the latest development run or of
# this PR's previous run.  On success, this run's database is in
# $RUNNER_TEMP/clang-tidy-db/out for the upload-artifact step.  Needs GH_TOKEN.
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

# Download the newest artifact NAME of a run on branch $1 of repository id $2 into in/$3.
fetch () {
    local run="" page resp
    # Newest first; looking at the newest 500 is plenty.
    for page in 1 2 3 4 5; do
        resp=$(gh api "repos/${GITHUB_REPOSITORY}/actions/artifacts?name=${NAME}&per_page=100&page=${page}") \
            || { echo "clang-tidy database: cannot list the artifacts ($3)"; return; }
        run=$(jq -r ".artifacts[] | select(.expired | not)
                     | select(.workflow_run.head_branch == \"$1\" and .workflow_run.head_repository_id == $2)
                     | .workflow_run.id" <<< "${resp}" | head -n 1)
        [[ -n "${run}" ]] && break
        [[ $(jq '.artifacts | length' <<< "${resp}") == 100 ]] || break
    done
    if [[ -z "${run}" ]]; then
        echo "clang-tidy database: no ${NAME} of branch $1 ($3)"
    elif gh run download "${run}" -R "${GITHUB_REPOSITORY}" -n "${NAME}" -D "${DB}/in/$3"; then
        echo "clang-tidy database: $3 run ${run}"
    else
        echo "clang-tidy database: cannot download run ${run} ($3)"
    fi
}

fetch development "${GITHUB_REPOSITORY_ID}" development
if [[ "${GITHUB_EVENT_NAME}" == "pull_request" ]]; then
    fetch "${GITHUB_HEAD_REF}" "$(jq .pull_request.head.repo.id "${GITHUB_EVENT_PATH}")" pr
fi

"${GITHUB_WORKSPACE}/Tools/C_scripts/mmclt.py" --input "${LOG}" \
    --db-in "${DB}"/in/*/*.txt --db-out "${DB}/out/${NAME}.txt" --db-dir "${DB}/hashes" \
    --clang-tidy="${CT}" --clang-tidy-args="${CT_ARGS}" || exit 1

make -j4 -k -f clang-tidy-ccache-misses.mak
