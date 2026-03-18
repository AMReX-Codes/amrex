# AMReX Agents Guide

Use this guide whenever you orchestrate explorers/workers inside the AMReX repository. It covers both AMReX developers (PR reviews, bug hunts, new features, and documentation) and AMReX users who ask agents for help learning or building with AMReX. AMReX itself is a C++/Fortran framework for block-structured adaptive mesh refinement (AMR) targeting large-scale PDE simulations on CPU and GPU architectures (CUDA, HIP, SYCL).

## Purpose & Personas

- **AMReX developers** – structure every agent task (reviews, fixes, features, documentation) so it is scoped, reproducible, and merged with confidence.
- **AMReX users** – route questions about capabilities, docs, tutorials, builds, or troubleshooting through the authoritative resources already shipped with the repo.

**Navigation tip**: We reference repository docs by section titles rather than line numbers. Use `rg -n "<heading text>" <file>` (or your editor’s outline) to jump to the relevant sections quickly.

## Repository Layout at a Glance

- `Src/` – Primary C++/Fortran implementation. Subfolders group functionality:
  - `Amr` / `AmrCore` house mesh hierarchy management, tagging, and regridding logic.
  - `Base` collects runtime essentials (memory arenas, geometry, I/O helpers) shared by every backend.
  - `Boundary`, `EB`, `LinearSolvers`, `Particle`, `FFT`, etc. provide focused subsystems; check their `CMakeLists.txt` for build toggles before touching code.
  - `Extern` and `F_Interfaces` bridge external packages and Fortran bindings.
- `Tests/` – CTest targets and sample drivers organized by topic (EB, GPU, LinearSolvers, Particles, etc.). When enabling `AMReX_ENABLE_TESTS`, these turn into runnable executables (`ctest -N` to list).
- `Docs/` – Source of all published documentation. The `sphinx_documentation/` tree feeds the public HTML docs; `Doxygen/` supports reference builds. Edit these when updating guides.
- `Tools/` – Build helpers, scripts, and shared CMake modules (e.g., `Tools/CMake/AMReXOptions.cmake`) and GNU makefiles.

## Operating Principles

- **Branch hygiene**: Work from short-lived branches based on the latest `development`, and never commit directly on the tracking `development` branch (see the “Git workflow” section of `CONTRIBUTING.md`, especially the “Generally speaking” rules about keeping `development` clean).
- **Single integration branch**: Treat `development` as the one authoritative branch AMReX maintains (see “Development Model” in `CONTRIBUTING.md`). Every PR must target it, monthly releases are tagged from it, and local work should always rebase onto it before review.
- **Coding style**: Follow the “AMReX Coding Style Guide” for indentation, brace usage, spacing, and member naming (see the “AMReX Coding Style Guide” section in `CONTRIBUTING.md`). If a change touches code and documentation, keep the style fixes local to the edited blocks.
- **Plan, scope, and delegate**: For any non-trivial task, sketch a plan, assign clear ownership when spawning explorers/workers, and avoid overlapping write scopes. Prefer `rg` for repo searches to stay fast in large trees.
- **Build and test defaults**: Confirm which build system the target supports. Repository-level libraries and executables use `cmake` as described in `Docs/sphinx_documentation/source/BuildingAMReX.rst` (section “Customization options”), with `ctest` as the default verification step and flags like `-DAMReX_ENABLE_TESTS=ON` plus `-DAMReX_TEST_TYPE=Small` when you only need a light signal (see the “Tests” block inside `Tools/CMake/AMReXOptions.cmake`). Most tests and tutorials also ship a `GNUmakefile`, and a few legacy drivers only expose that path, so `cd` into the test directory and run `make -j` with the variables it expects (e.g., `DIM`, `USE_MPI`, `USE_CUDA`, `COMP`) following `Tools/GNUMake/README.md`.
- **Documentation sources**: Lean on the curated entry points listed in the “Documentation” section of `README.md`. They point to the public Sphinx build at `https://amrex-codes.github.io/amrex/docs_html/`, which mirrors the sources under `Docs/sphinx_documentation`. Treat the standalone tutorials repository referenced in `Tutorials/README.md` (`https://github.com/AMReX-Codes/amrex-tutorials`) as additional runnable examples you can cite when users ask how to get started.
- **Issue logging & hand-off**: Keep a personal, untracked scratchpad on each machine (we recommend `agent-notes/<NN>-<component>-<short-description>.md`). Use it to capture open questions, repro notes, or follow-ups, reusing the numbering/component/title convention described below. Include suggested patches whenever possible so the next agent can act quickly.
- **Learn from past bugs**: If you already keep a local `agent-notes/` notebook, skim it before diving into similar code to refresh common pitfalls—many historical AMReX bugs came from copy-paste mistakes (e.g., duplicated kernels, swapped indices, missing constant updates), so assume near-identical blocks may hide divergences.

## Developer Playbooks

### PR and Bug Reviews
1. **Sync & inspect** – Update the local branch, note the PR/issue scope, and record file ownership expectations.
2. **Reproduce & read** – Reproduce the report using the author’s steps or by running the focused tests. While reading diffs, confirm they honor the rules in the “AMReX Coding Style Guide” section of `CONTRIBUTING.md`.
3. **Hunt for copy-paste drift** – Compare mirrored kernels, dimension-specific code paths, and duplicated tables; historical regressions often stem from edits applied to one block but not its sibling. Look for suspiciously similar snippets that differ only by variable names or miss a constant update.
4. **Verify** – Configure the project with the appropriate options (for example, GPU flags or `-DAMReX_TEST_TYPE=Small`) and run `ctest --output-on-failure` from the build directory, or use the test’s `make` rule when it only ships a `GNUmakefile`.
5. **Focus on hot spots** – Use `cmake --build build -j --target <target_name>` for a single executable/test and `ctest --test-dir build -R <regex>` (or `ctest -R <regex>` inside the build tree) to rerun only the impacted cases; for `GNUmakefile` flows, rerun `make -j` (optionally with a target such as `make run` or `make tests`) inside the test directory. Capture the output.
6. **Report** – Summarize findings (blocking issues first), highlight required tests, and cite files/lines that need attention.
7. **Log follow-ups** – If more work is required, open or update the matching file in `agent-notes/` (or your local scratchpad) so the next agent inherits context.

### Feature or Fix Implementation
1. **Understand scope** – Capture requirements, physics context, and success criteria from the originating issue/PR.
2. **Configure builds quickly** – Choose the workflow the directory expects. Use the standard `cmake` pattern below, adding any extra `-D` knobs listed in the “Customization options” portion of `Docs/sphinx_documentation/source/BuildingAMReX.rst`.

   ```bash
   cmake -S . -B build \
     -DAMReX_ENABLE_TESTS=ON \
     -DAMReX_TEST_TYPE=Small
   cmake --build build -j
   ctest --test-dir build --output-on-failure
   ```

   When only one binary or test matters, leverage `cmake --build build -j --target <target_name>` and `ctest --test-dir build -R <regex>` to keep feedback loops short.

   Directories that rely on `GNUmakefile` (many tutorials/tests, plus a handful of legacy drivers) follow the guidance in `Docs/sphinx_documentation/source/BuildingAMReX.rst` and `Tools/GNUMake/README.md`. Set only the variables that the specific example requires (`DIM`, `USE_MPI`, `USE_CUDA`, `COMP`, etc.) so they reflect the hardware/features you intend to exercise. Edit the local `GNUmakefile` or pass those variables on the command line, then build with `make`. For instance, a 3D CNS run that enables both MPI and CUDA would be:

   ```bash
   cd Tests/GPU/CNS
   make -j8 DIM=3 USE_MPI=TRUE USE_CUDA=TRUE
   ```

3. **Implement with traceability** – Touch only the files you own in this task, annotate complex code with succinct comments, and reference relevant issue IDs.
4. **Document** – Update user-facing docs whenever behavior changes. Pull content from the “Documentation” section of `README.md` (User’s Guide, Example Codes, Guided Tutorials, Technical Reference) so users know where to look.
5. **Hand off** – Record remaining questions, test logs, or benchmarking data inside `agent-notes/` (or your local scratchpad) or the PR description, including exact commands run and their outcomes.

### Documentation and Tutorial Updates
- For feature additions, mirror the doc hierarchy described in the “Documentation” section of `README.md` so the User’s Guide, Example Codes, and Guided Tutorials stay synchronized.
- Surface new build options or workflows in `Docs/sphinx_documentation/source/BuildingAMReX.rst` so the `Customization options` table stays authoritative.

## Guidance for AMReX Users Working with Agents

- **Getting oriented**: Summarize AMReX capabilities using the “Overview,” “Features,” and “Documentation” sections in `README.md`. Link users to the appropriate resource (User’s Guide, Example Codes, Guided Tutorials, Technical Reference).
- **Building & testing quickly**: Start with whichever build system the example ships. Walk users through the `cmake` workflow highlighted in the “Customization options” part of `Docs/sphinx_documentation/source/BuildingAMReX.rst` (showing how to toggle features with the `-D<var>=<value>` syntax), and when they are inside a tutorial or test directory that provides a `GNUmakefile`, point them to the same doc plus `Tools/GNUMake/README.md` so they can run `make -j` with variables such as `DIM`, `USE_MPI`, and `USE_CUDA`.
- **Learning resources**: Direct users to the standalone tutorials repository noted in `Tutorials/README.md` (`https://github.com/AMReX-Codes/amrex-tutorials`) and supplement with the slides/videos featured near the “Documentation” section of `README.md`.
- **Consult Sphinx sources**: When clarifying documentation or preparing local updates, read directly from `Docs/sphinx_documentation` (especially the `source/` subtree). This is the exact content published online, so citing it keeps agent answers aligned with the official docs.
- **Getting help or contributing back**: Encourage questions through GitHub Discussions and remind users that contributions go through `CONTRIBUTING.md`, as described in the “Get Help” and “Contribute” sections of `README.md`.

## Agent Notes & Hand-off

Agents rely on a lightweight, per-machine scratchpad to capture ephemeral context (repro steps, local experiments, or future TODOs) without polluting the repo. This is an agent-side convention, not an upstream AMReX requirement—keep it untracked so you can jot candid notes and prune freely.

- **Where**: Create an `agent-notes/` folder at the repo root. If you already standardized on another name, that’s acceptable—just stay consistent on that machine. File names follow `NN-component-short-description.md`, where `NN` is a zero-padded counter unique per workstation.
- **What to include**:
  - Title line summarizing the issue or follow-up.
  - Metadata bullets for `Type` (Bug/Feature/Docs), `Severity`, `Component`, and an approximate `Location` (file:line or directory).
  - Sections for `Problem`, `Impact`, and `Next steps` or `Suggested patch`. Link to relevant PRs, branches, or external tickets if applicable.
  - Exact reproduce/build/test commands and outputs to save the next agent time.
- **Sharing**: Because the folder is local-only, copy the relevant markdown snippet into a PR description, long-form review, or upstream issue whenever collaborators need visibility.

Include ready-to-apply patches or diff hunks whenever possible so other agents (or future you) can fast-track the fix.

## Quick Checklist

1. Confirm you are on a task-specific branch that tracks `development` cleanly (see the “Git workflow” guidance in `CONTRIBUTING.md`).
2. Plan the task, noting deliverables, ownership, and validation steps before spawning sub-agents.
3. Build with the workflow the directory expects: either run the standard `cmake`/`ctest` flow (with `AMReX_ENABLE_TESTS` and `AMReX_TEST_TYPE` toggles per “Customization options” in `Docs/sphinx_documentation/source/BuildingAMReX.rst` and the “Tests” block in `Tools/CMake/AMReXOptions.cmake`) or `cd` into the `GNUmakefile` tree and run `make -j` with the required variables (e.g., `DIM`, `USE_MPI`, `USE_CUDA`).
4. Update documentation and user guidance by referencing the resources enumerated in the “Documentation” section of `README.md`.
5. Capture unresolved work, context, and suggested patches in `agent-notes/` (or your local scratchpad) so future agents can pick up where you left off.
