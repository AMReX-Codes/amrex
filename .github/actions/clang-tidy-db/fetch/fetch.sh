#!/usr/bin/env bash
#
# Download this job's clang-tidy databases of the latest development run and
# of this PR's previous run into $RUNNER_TEMP/clang-tidy-db/in for
# clang-tidy-ccache-misses.sh.  Needs GH_TOKEN with actions: read.

set -u

NAME=clang-tidy-db-${GITHUB_WORKFLOW}-${GITHUB_JOB}
DB=${RUNNER_TEMP}/clang-tidy-db
mkdir -p "${DB}/in"

# Download the newest artifact NAME of a run on branch $1 of repository id $2
# into in/$3.  Normally two API calls: the list of artifacts named NAME,
# newest first, and the zip.
fetch () {
    local found="" page resp
    for page in 1 2 3; do  # 300 runs of this workflow is plenty
        resp=$(gh api "repos/${GITHUB_REPOSITORY}/actions/artifacts?name=${NAME}&per_page=100&page=${page}") \
            || { echo "clang-tidy database: cannot list the artifacts ($3)"; return; }
        found=$(jq -r ".artifacts[] | select(.expired | not)
                       | select(.workflow_run.head_branch == \"$1\" and .workflow_run.head_repository_id == $2)
                       | \"\(.id) \(.workflow_run.id)\"" <<< "${resp}" | head -n 1)
        [[ -n "${found}" ]] && break
        [[ $(jq '.artifacts | length' <<< "${resp}") == 100 ]] || break
    done
    if [[ -z "${found}" ]]; then
        echo "clang-tidy database: no ${NAME} of branch $1 ($3)"
        return
    fi
    local id=${found% *} run=${found#* }
    if gh api "repos/${GITHUB_REPOSITORY}/actions/artifacts/${id}/zip" > "${DB}/$3.zip" \
        && unzip -q -o "${DB}/$3.zip" -d "${DB}/in/$3"; then
        echo "clang-tidy database: $3 run ${run}"
    else
        echo "clang-tidy database: cannot download run ${run} ($3)"
    fi
}

fetch development "${GITHUB_REPOSITORY_ID}" development
if [[ "${GITHUB_EVENT_NAME}" == "pull_request" ]]; then
    fetch "${GITHUB_HEAD_REF}" "$(jq .pull_request.head.repo.id "${GITHUB_EVENT_PATH}")" pr
fi
