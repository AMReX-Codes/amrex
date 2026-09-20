# Security Policy

## Reporting security issues

To report a security issue, please follow the GitHub instructions on
privately reporting a security vulnerability
(https://docs.github.com/en/code-security/security-advisories/guidance-on-reporting-and-writing-information-about-vulnerabilities/privately-reporting-a-security-vulnerability#privately-reporting-a-security-vulnerability).

Please report privately rather than opening a public issue, so that a fix can
be prepared before the problem is widely known. A report is most useful if it
says which AMReX version is affected, how to reproduce the problem, and what
impact you believe it has.

For context on what AMReX treats as a security issue, and where its trust
boundaries lie, see the assurance case below. Because
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

## Security assurance case

This section sets out why AMReX's security requirements are met: what AMReX is
exposed to, where its trust boundaries lie, which design principles apply, and
how common implementation weaknesses are countered.

### Scope and threat model

AMReX is a library. It is compiled into an application that the user builds and
runs themselves, normally on a workstation or on an HPC allocation they already
hold. It is not a service, has no daemon or long-running process of its own,
does not listen on a network socket, does not authenticate users, and stores no
credentials.

The inputs AMReX processes -- inputs files, runtime parameters, plotfiles and
checkpoint files -- are supplied by the same person running the simulation, and
are read with the privileges that person already has. AMReX therefore does not
treat its own inputs as hostile: an input able to make AMReX misbehave grants
the author nothing they could not do directly by running their own code.

The realistic threats are consequently:

1. **Supply chain.** An attacker modifies AMReX source, a release artifact, or a
   build-time dependency so that applications built from it are compromised.
2. **Memory-safety defects.** A bug in AMReX causes out-of-bounds access or
   undefined behavior, which is primarily a correctness and reliability problem
   but could be reachable from data an application accepts from elsewhere.
3. **Downstream misuse.** An application built on AMReX exposes AMReX to data
   from a less trusted source than AMReX itself assumes.

Threat 3 is outside AMReX's control; applications that read untrusted files are
responsible for validating them before handing them to AMReX. This document
addresses threats 1 and 2.

### Trust boundaries

AMReX crosses no privilege boundary at runtime. Code and data enter the same
trust domain they came from:

| Interface | Trust |
|---|---|
| Application code calling AMReX | Same process, same privileges -- no boundary |
| Inputs files, runtime parameters | Supplied by the user running the job |
| Plotfile and checkpoint I/O | Written and read by the same user |
| MPI communication | Within one job on a trusted interconnect |
| GPU kernels (CUDA, HIP, SYCL) | Same user, same allocation |

The boundary that does matter is at development and distribution time: the
GitHub repository, the CI that builds and tests it, and the release artifacts
users download.

### Secure design principles

- **Least privilege.** AMReX requires no elevated rights. Nothing in it is
  installed setuid, runs as a service, or asks for privileges beyond those of
  the invoking user.
- **Small attack surface.** No network listeners, no credential handling, no
  cryptography, no plugin loading from untrusted paths.
- **Economy of mechanism.** Optional dependencies (MPI, HDF5, HYPRE, PETSc,
  SUNDIALS and others) are found at build time with CMake `find_package` and
  supplied by the system or by Spack/conda. AMReX vendors no forked copies of
  external libraries, so a security update to any of them takes effect without
  changes to AMReX.
- **Defense of the supply chain.** Changes reach `development` through reviewed
  pull requests, CI runs on every one of them, and Dependabot tracks the GitHub
  Actions the project depends on.

### Countering common implementation weaknesses

The OWASP Top 10 is largely inapplicable: AMReX has no web surface, no
authentication, no session handling, and no injection sinks of the kind it
describes. The weaknesses that do apply are the memory-safety and
undefined-behavior classes of CWE, countered as follows.

| Measure | Where |
|---|---|
| Static analysis on every push and pull request | CodeQL, `.github/workflows/codeql.yml` |
| clang-tidy, including `modernize-*` and `deprecated-declarations` | `.clang-tidy`, run in seven CI workflows |
| Warnings treated as errors | `-Werror` in `Tools/CMake/AMReXFlagsTargets.cmake` and `Tools/GNUMake/comps/gnu.mak` |
| Address, undefined-behavior and thread sanitizer builds | `FSANITIZER=TRUE` and `THREAD_SANITIZER=TRUE` in the GNU Make build |
| Assertions, on by default in debug builds | `AMReX_ASSERTIONS`, `AMREX_ASSERT` / `AMREX_ALWAYS_ASSERT` |
| Optional bounds checking on array access | `AMReX_BOUND_CHECK` for `Array4` |
| Continuous build and test across compilers and GPU backends | GCC, Clang, Intel, CUDA, HIP, SYCL, macOS, Windows |
| Periodic deep static analysis | Coverity Scan |
| Modern C++ | C++20 required; raw owning pointers and C arrays discouraged by the clang-tidy configuration |

Memory safety is not guaranteed -- AMReX is C++ and performance-critical, and
bounds checking and sanitizers are opt-in rather than always on. The argument
here is that defects of this class are systematically hunted rather than
eliminated by construction, and that their impact is bounded by the trust model
above: AMReX runs on its own user's data with its own user's privileges.
