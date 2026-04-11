# AMReX Agents Guide

Use this guide whenever you orchestrate explorers/workers inside the AMReX repository. It keeps PR reviews, bug hunts, new features, docs, and user support consistent for contributors and AI copilots alike.

## Purpose & Personas

- **AMReX developers** – follow a predictable workflow so fixes/features land quickly and safely.
- **AMReX users** – receive accurate build/test/docs guidance sourced from the repo.

## Repository Map

- `Src/` – Core C++/Fortran implementation. `Amr*` manage hierarchy/regridding, `Base` hosts runtime utilities, and folders such as `Boundary`, `EB`, `LinearSolvers`, `Particle`, `FFT` contain subsystem code (check local `CMakeLists.txt` for toggles).
- `Tests/` – Organized by topic; when `AMReX_ENABLE_TESTS=ON`, they expose `ctest` targets. Many tutorials/tests also ship a `GNUmakefile` for standalone runs.
- `Docs/` – Authoritative Sphinx sources under `Docs/sphinx_documentation/`; Doxygen lives under `Docs/Doxygen/`.
- `Tools/` – Shared CMake modules, GNU make helpers, scripts.

## Build & Test Reference

- **Default CMake flow**:
  ```bash
  cmake -S . -B build \
    -DAMReX_ENABLE_TESTS=ON \
    -DAMReX_TEST_TYPE=Small
  cmake --build build -j
  ctest --test-dir build --output-on-failure
  ```
  Toggle options from `Docs/sphinx_documentation/source/BuildingAMReX.rst` (“Customization options”) and `Tools/CMake/AMReXOptions.cmake` for GPUs, dimensions, debug flags, etc.
- **Targeted builds/tests**: `cmake --build build -j --target <name>` for individual executables; `ctest --test-dir build -R <regex>` (or `ctest -R <regex>`) to rerun impacted cases only.
- **GNUmakefile workflows**: When a directory ships a `GNUmakefile`, `cd` there and run `make -j` with required variables (e.g., `DIM`, `USE_MPI`, `USE_CUDA`, `COMP`) as documented in `Docs/sphinx_documentation/source/BuildingAMReX.rst` and `Tools/GNUMake/README.md`.
- **Log everything**: Capture exact commands plus pass/fail output in PR descriptions or linked issues so reviewers can reproduce without guessing.

## Hard Rules & Defaults

- Work on short-lived topic branches; never commit directly to `development`. Rebase before opening or updating a PR.
- Plan before you edit. Define scope, owners, and validation steps, and use `rg` (or your editor) for fast searches to avoid wandering.
- Choose a build workflow from the “Build & Test Reference” and stick to it for the task; don’t restate commands elsewhere—link back and record the actual invocations you ran.
- Update user-facing docs under `Docs/sphinx_documentation/` whenever behavior changes, and mention the doc edits in the same PR.
- When adding Doxygen comments or Sphinx docs, write for AMReX users: explain behavior and usage, omit unnecessary implementation details, and keep the tone instructional instead of developer-oriented.
- Never log work inside `CHANGES.md`; release notes are curated separately, so leave that file untouched unless maintainers ask for a release update.
- Hand-off context inside PRs/issues (commands run, results, TODOs). Personal scratchpads are optional, local, and should be pruned frequently.
- AI delegation: assign disjoint write scopes, keep subtasks bounded, and serialize conflicting edits when necessary. Review every AI-generated diff like a teammate’s patch.
- Safety: never paste secrets, private datasets, or PII into prompts; scrub logs before sharing; cite sources and licensing when in doubt by escalating to maintainers.
- Testing accountability: every substantive change must state which tests ran (and command lines) in the PR/issue.
- Historical hotspots: mirrored kernels and dimension-specific paths often drift—compare siblings whenever you touch them.

## Task Playbooks

### PR & Bug Review
1. **Sync & scope** – Update the branch, read the PR/issue description, and note ownership expectations.
2. **Reproduce** – Follow the reporter’s steps or minimal `ctest`/`make` invocations from the Build & Test Reference.
3. **Inspect diffs** – Check style, compare mirrored kernels/dimensional variants, and look for missing constant updates.
4. **Verify** – Run the focused tests (GPU flags, `-DAMReX_TEST_TYPE=Small`, GNU make targets, etc.) and capture output.
5. **Report** – List blockers first, cite files/lines, and record the commands/tests in the PR.
6. **Follow-ups** – Document remaining work in the PR/issue; use a scratchpad only for personal reminders.

### Feature or Fix Implementation
1. **Confirm requirements** – Capture physics context, success metrics, and acceptance tests from the originating ticket.
2. **Configure builds** – Apply the Build & Test Reference workflow that matches the directory; use targeted targets/tests for fast iteration.
3. **Implement carefully** – Touch only owned files, note tricky code with short comments, and maintain coding-style guidelines.
4. **Update docs** – When behavior or UX changes, add/adjust content in the relevant `Docs/sphinx_documentation/` section (User’s Guide, Technical Reference, tutorials, etc.).
5. **Validate** – Run and record the necessary `ctest`/`make` commands plus any benchmarking or profiling relevant to the change.
6. **Hand off** – Summarize status, remaining risks, and next steps in the PR description; attach log excerpts only if scrubbed.

### Documentation & Tutorial Updates
- Mirror the hierarchy noted in `README.md` (User’s Guide, Example Codes, Guided Tutorials, Technical Reference) so published docs stay synchronized.
- Surface any new build knobs/workflows in `Docs/sphinx_documentation/source/BuildingAMReX.rst`.
- Reference runnable examples in `Tutorials/README.md` or `https://github.com/AMReX-Codes/amrex-tutorials` when guiding users.

## Guidance for AMReX Users Working with Agents

- **Orientation** – Summarize capabilities using the “Overview,” “Features,” and “Documentation” sections of `README.md`, then link users to the best resource (User’s Guide, Example Codes, Guided Tutorials, Technical Reference).
- **Build help** – Walk through the Build & Test Reference commands, clarifying how to toggle features via `-DVAR=value` or `make` variables the tutorial expects.
- **Learning resources** – Point to the tutorials repo plus any slides/videos referenced near the Documentation section of `README.md`.
- **Authoritative sources** – Read directly from `Docs/sphinx_documentation/` when answering doc questions so citations match the published site.
- **Support channels** – Encourage GitHub Discussions/issues for unresolved questions and remind users that contributions follow `CONTRIBUTING.md`.

## Optional Personal Notes

A local scratchpad can help you capture quick reminders, but keep it untracked, short-lived, and private. Anything others need to know belongs in PR descriptions, issues, or review comments.

## Quick Checklist

1. On a topic branch rebased on `development`? If not, fix it.
2. Do you have a written plan that names owners, validation, and Build & Test commands? If not, write one before editing.
3. Are tests/docs updated and the exact commands/results logged in the PR/issue? If not, add them.
4. Are delegation, safety, and hand-off notes captured in canonical threads (not just scratchpads)? If not, update them now.
