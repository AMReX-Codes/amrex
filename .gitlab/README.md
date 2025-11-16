This note describes the setup of our GitHub/GitLab integration. The main
AMReX repo is hosted at https://github.com/AMReX-Codes/amrex, with a mirror
at https://gitlab.spack.io/amrex/amrex. The GitLab mirror runs CI tests on
runners hosted at
https://systems.nic.uoregon.edu/internal-wiki/index.php?title=Category:Servers. We
run three types of tests: nightly tests using the development branch, tests
triggerd when new commits are pushed to the GitLab mirror (either through
GitLab's regular pulling of the latest development branch from GitHub or by
the developers directly), and tests triggered by a special comment from
maintainers.

Instructions for creating a pull mirror can be found at
https://docs.gitlab.com/user/project/repository/mirror/pull/. You can choose
to pull selected branches or all branches. For AMReX, this choice does not
matter since we only have one active branch. However, for repoositories with
many active branches, pulling all branches may trigger CI jobs on every
activity.

Instructions for creating schedules pipelines can be found at
https://docs.gitlab.com/ci/pipelines/schedules/. By default, GitLab uses
`.gitlab-ci.yml` at the repository root for pipeline configurations, but
this can be changed. AMReX's GitLab CI file is in the `.gitlab`
directory. To change the configuration file location, select `Settings ->
CI/CD -> General pipelines` and update `CI/CD configuration file`.

CI/CD jobs triggered by push require no special setup. They can also be
disabled.

For GitHub PR triggered jobs, we choose not to run automatically. Instead
they need to be triggered by a `/run-hpsf-gitlab-ci` comment from authorized
maintainers, which starts the `run-hpsf-gitlab-ci` workflow defined in
`.github/workflows/trigger-hpsf-gitlab-ci.yml`. This workflow uses GitLab's
REST API to start a pipeline job on GitLab. The PR comment triggered
pipeline job uses the same configuration file shared with the scheduled
pipeline: `.gitlab/hpsf-gitlab-ci.yml`. The PR comment triggered job pulls
the PR branch from GitHub first before running tests. For this approach to
work, we store a pipeline trigger (obtained from GitLab's `Settings -> CI/CD
-> Pipeline trigger tokens`) as a secret at GitHub's `Settings -> Secrets
and variables -> Actions -> Repository secrets`. The GitHub workflow waits
for the result of the GitLab pipeline result and posts the final status and
a link to the result as a comment.

There is one complication with this approach. The GitHub PR triggered
pipeline appears on GitLab under the title of the latest commit of the
default branch (i.e., development), which can be confusing. To address this,
we add a step in the GitHub workflow to change the GitLab pipeline's name to
the PR's title. This step requires an API token, which diffs from the
pipeline trigger token. A project API token can be created from the
project's `Settings -> Access tokens`. It needs `api` scopes and the
`Maintainer` role. The token will expire in one year. After creating a key,
you must copy and store it on GitHub before leaving the page, as it will
become invisible afterward.
