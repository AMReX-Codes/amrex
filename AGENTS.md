# AMReX Agents Guide

Use this guide whenever you orchestrate explorers/workers inside the AMReX repository. It covers both AMReX developers (PR reviews, bug hunts, new features, and documentation) and AMReX users who ask agents for help learning or building with AMReX.

## Purpose & Personas

- **AMReX developers** – structure every agent task (reviews, fixes, features, documentation) so it is scoped, reproducible, and merged with confidence.
- **AMReX users** – route questions about capabilities, docs, tutorials, builds, or troubleshooting through the authoritative resources already shipped with the repo.

## Operating Principles

- **Branch hygiene**: Work from short-lived branches based on the latest `development`, and never commit directly on the tracking `development` branch (see `CONTRIBUTING.md:191-208`).
- **Single integration branch**: Treat `development` as the one authoritative branch AMReX maintains (see `CONTRIBUTING.md:1-35`). Every PR must target it, monthly releases are tagged from it, and local work should always rebase onto it before review.
- **Coding style**: Follow the “AMReX Coding Style Guide” for indentation, brace usage, spacing, and member naming (see `CONTRIBUTING.md:210-255`). If a change touches code and documentation, keep the style fixes local to the edited blocks.
- **Plan, scope, and delegate**: For any non-trivial task, sketch a plan, assign clear ownership when spawning explorers/workers, and avoid overlapping write scopes. Prefer `rg` for repo searches to stay fast in large trees.
- **Build and test defaults**: Configure builds with `cmake` as described in `Docs/sphinx_documentation/source/BuildingAMReX.rst:381-416`, and treat `CTest` as the default verification step. Enable tests via `-DAMReX_ENABLE_TESTS=ON` and use `-DAMReX_TEST_TYPE=Small` when a light signal is enough (`Tools/CMake/AMReXOptions.cmake:515-526`).
- **Documentation sources**: Lean on the curated entry points listed in `README.md:67-99` whenever users ask for learning material, examples, or contribution guidance, and link to tutorials mentioned in `Tutorials/README.md:1`.
- **Issue logging & hand-off**: Keep a personal (untracked) `issues/<NN>-<component>-<short-description>.md` directory on each machine. Use it to capture open questions, repro notes, or follow-ups, reusing the numbering/component/title convention described below. Include suggested patches whenever possible so the next agent can act quickly.
- **Learn from past bugs**: Before diving into similar code, skim your local `issues/` notes to refresh common pitfalls—many historical AMReX bugs came from copy-paste mistakes (e.g., duplicated kernels, swapped indices, missing constant updates), so assume near-identical blocks may hide divergences.

## Developer Playbooks

### PR and Bug Reviews
1. **Sync & inspect** – Update the local branch, note the PR/issue scope, and record file ownership expectations.
2. **Reproduce & read** – Reproduce the report using the author’s steps or by running the focused tests. While reading diffs, confirm they honor the style guide (`CONTRIBUTING.md:210-255`).
3. **Hunt for copy-paste drift** – Compare mirrored kernels, dimension-specific code paths, and duplicated tables; historical regressions often stem from edits applied to one block but not its sibling. Look for suspiciously similar snippets that differ only by variable names or miss a constant update.
4. **Verify** – Configure the project with the appropriate options (for example, GPU flags or `-DAMReX_TEST_TYPE=Small`) and run `ctest --output-on-failure` from the build directory.
5. **Focus on hot spots** – Use `cmake --build build -j --target <target_name>` for a single executable/test and `ctest --test-dir build -R <regex>` (or `ctest -R <regex>` inside the build tree) to rerun only the impacted cases, then capture the output.
6. **Report** – Summarize findings (blocking issues first), highlight required tests, and cite files/lines that need attention.
7. **Log follow-ups** – If more work is required, open or update the matching file under `issues/` so the next agent inherits context.

### Feature or Fix Implementation
1. **Understand scope** – Capture requirements, physics context, and success criteria from the originating issue/PR.
2. **Configure builds quickly** – Use the standard pattern below, adding any extra `-D` knobs listed in `Docs/sphinx_documentation/source/BuildingAMReX.rst:381-446`.

   ```bash
   cmake -S . -B build \
     -DAMReX_ENABLE_TESTS=ON \
     -DAMReX_TEST_TYPE=Small
   cmake --build build -j
   ctest --test-dir build --output-on-failure
   ```

   When only one binary or test matters, leverage `cmake --build build -j --target <target_name>` and `ctest --test-dir build -R <regex>` to keep feedback loops short.

3. **Implement with traceability** – Touch only the files you own in this task, annotate complex code with succinct comments, and reference relevant issue IDs.
4. **Document** – Update user-facing docs whenever behavior changes. Pull content from the `Documentation` section of `README.md:67-74` (User’s Guide, Example Codes, Guided Tutorials, Technical Reference) so users know where to look.
5. **Hand off** – Record remaining questions, test logs, or benchmarking data inside `issues/` or the PR description, including exact commands run and their outcomes.

### Documentation and Tutorial Updates
- For feature additions, mirror the doc hierarchy described in `README.md:67-99` so the User’s Guide, Example Codes, and Guided Tutorials stay synchronized.
- When tutorials move or expand, note the canonical source (`Tutorials/README.md:1` points to the external tutorials repository) and update both AMReX docs and the tutorial repo as needed.
- Surface new build options or workflows in `Docs/sphinx_documentation/source/BuildingAMReX.rst` so the `Customization options` table stays authoritative.

## Guidance for AMReX Users Working with Agents

- **Getting oriented**: Summarize AMReX capabilities using the “Overview,” “Features,” and “Documentation” sections in `README.md:51-99`. Link users to the appropriate resource (User’s Guide, Example Codes, Guided Tutorials, Technical Reference).
- **Building & testing quickly**: Walk users through the `cmake` workflow highlighted in `Docs/sphinx_documentation/source/BuildingAMReX.rst:381-416`, and explain how to toggle features with the `-D<var>=<value>` syntax using examples from the same section.
- **Learning resources**: Point to tutorials via `Tutorials/README.md:1`, plus any relevant slides or videos referenced in `README.md:69-74`.
- **Consult Sphinx sources**: When clarifying documentation or preparing local updates, read directly from `Docs/sphinx_documentation` (especially the `source/` subtree). This is the exact content published online, so citing it keeps agent answers aligned with the official docs.
- **Getting help or contributing back**: Encourage questions through GitHub Discussions and remind users that contributions go through `CONTRIBUTING.md` as noted in `README.md:87-99`.

## Issues & Hand-off Notes

- **Where**: Create an `issues/` folder at the repo root on every workstation (it remains untracked). File names follow `NN-component-short-description.md`, where `NN` is a zero-padded counter that stays unique per machine.
- **What to include**:
  - Title line summarizing the issue or follow-up.
  - Metadata bullets for `Type` (Bug/Feature/Docs), `Severity`, `Component`, and an approximate `Location` (file:line or directory).
  - Sections for `Problem`, `Impact`, and `Next steps` or `Suggested patch`. Link to relevant PRs, branches, or external tickets if applicable.
  - Exact reproduce/build/test commands and outputs to save the next agent time.
- **Sharing**: Because the folder is local-only, copy the relevant markdown snippet into a PR description, long-form review, or upstream issue whenever collaborators need visibility.

Include ready-to-apply patches or diff hunks whenever possible so other agents (or future you) can fast-track the fix.

## Quick Checklist

1. Confirm you are on a task-specific branch that tracks `development` cleanly (`CONTRIBUTING.md:191-208`).
2. Plan the task, noting deliverables, ownership, and validation steps before spawning sub-agents.
3. Build with the right `cmake` options and run `ctest` (using `AMReX_ENABLE_TESTS` and `AMReX_TEST_TYPE` toggles) to validate changes (`Docs/sphinx_documentation/source/BuildingAMReX.rst:381-416`, `Tools/CMake/AMReXOptions.cmake:515-526`).
4. Update documentation and user guidance by referencing the resources enumerated in `README.md:67-99` and `Tutorials/README.md:1`.
5. Capture unresolved work, context, and suggested patches in `issues/` so future agents can pick up where you left off.
