# AMReX Security Assurance Case

This document sets out why AMReX's security requirements are met: what AMReX is
exposed to, where its trust boundaries lie, which design principles apply, and
how common implementation weaknesses are countered. To report a vulnerability,
see [SECURITY.md](SECURITY.md).

## Scope and threat model

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

## Trust boundaries

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

## Secure design principles

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

## Countering common implementation weaknesses

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

## Reporting and response

Vulnerabilities are reported privately through GitHub's security advisory
mechanism, as described in [SECURITY.md](SECURITY.md). The technical committee
listed in [GOVERNANCE.md](GOVERNANCE.md) is responsible for triage and fixes.
