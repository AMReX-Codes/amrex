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
and variables -> Actions -> Repository secrets`.

After the GitLab pipeline finishes, its `.post` stage will post the final
status and a link to the results back to the GitHub PR as a comment. This is
done through a GitHub App that we built and installed in the AMReX
repository. The App was created via
https://github.com/organizations/AMReX-Codes/settings/apps. It does not need
"Webhook" access. For repository permissions, it only needs read & write to
pull requests so it can create PR comments. The app requires a private key,
which you generate during the setup. After the app was installed in the
amrex repository, we got an installation ID. We then stored the app ID,
installation ID and the private key in GitLab's `Settings -> CI/CD ->
Variables`. The app ID isn't a secret. So you can store it as clear text. In
fact, GitLab does not allow 7-digit masked variables anyway. The
installation ID is also not senstive, but nevertheless we stored it as
protected and masked. The private key is a secret that must be protected and
masked. We also diabled "Expand" for all of these variables because the CI
script doesn't need variable expansion. GitLab seems to have a bug that
prevents saving the private key as a multi-line value, and saving it as a
file didn't work either. So we encoded it with `base64 -w0` to turn it into
a single line. That's why in the GitLab CI script we have to decode it.
