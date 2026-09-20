# Reporting Security Issues

To report a security issue, please follow the GitHub instructions on
privately reporting a security vulnerability
(https://docs.github.com/en/code-security/security-advisories/guidance-on-reporting-and-writing-information-about-vulnerabilities/privately-reporting-a-security-vulnerability#privately-reporting-a-security-vulnerability).

Please report privately rather than opening a public issue, so that a fix can
be prepared before the problem is widely known. A report is most useful if it
says which AMReX version is affected, how to reproduce the problem, and what
impact you believe it has.

For context on what AMReX treats as a security issue, and where its trust
boundaries lie, see [SECURITY_ASSURANCE.md](SECURITY_ASSURANCE.md). Because
AMReX is a library that runs with the privileges of the user who launched it,
and reads inputs supplied by that same user, many defects that would be
vulnerabilities in a networked service are ordinary bugs here. Report anything
you are unsure about privately and we will judge it together; correctness bugs
can be reported normally through GitHub Issues.

## How we respond

Reports are received by the AMReX technical committee, whose members and
responsibilities are listed in [GOVERNANCE.md](GOVERNANCE.md).

1. **Acknowledgement.** We aim to acknowledge a report within a few business
   days, and to tell you whether we consider it a security issue or an ordinary
   bug.
2. **Assessment.** We reproduce the problem where we can, determine which
   versions are affected, and agree on severity with the reporter.
3. **Fix.** A fix is prepared, reviewed, and merged into `development`. AMReX
   tags a release on the first workday of each month; a fix judged serious
   enough may prompt a release outside that schedule.
4. **Disclosure.** Once a fix is available we publish a GitHub security
   advisory describing the problem, the affected versions and the fix.
5. **Credit.** Reporters are credited in the advisory unless they ask not to
   be.

If a report turns out not to be a security issue, we will say so and, with your
agreement, continue handling it as a normal bug report in public.
