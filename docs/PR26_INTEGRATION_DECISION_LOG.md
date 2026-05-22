# PR 26 Integration Decision Log

This document summarizes the PR #26 integration work and decisions made so far.
It is intended as a compact handover for future Codex sessions. Use it together
with the detailed task checklist in
[`docs/PR26_INTEGRATION_HANDOVER_PLAN.md`](./PR26_INTEGRATION_HANDOVER_PLAN.md).

## Current State

- PR being integrated: GitHub PR #26 for `PDBeurope/AFDB-Integration-Kit`.
- PR snapshot branch: `review/pr-26`.
- Parent integration branch: `integration-pr-26-gpu`.
- Current active branch at the end of this conversation:
  `integration-pr-26-gpu-step-4-dssp`.
- Step branches completed so far:
  - `integration-pr-26-gpu-step-1-dependencies`
  - `integration-pr-26-gpu-step-2-hygiene`
  - `integration-pr-26-gpu-step-3-install-docs`
  - `integration-pr-26-gpu-step-4-dssp`
- Merged into `integration-pr-26-gpu` so far:
  - Step 1 dependencies/test baseline
  - Step 2 mechanical hygiene
  - Step 3 install documentation
- Not yet merged into `integration-pr-26-gpu` at the end of this conversation:
  - Step 4 DSSP refactor branch

The immediate next action should be to merge Step 4 back into
`integration-pr-26-gpu`, verify the parent branch, and then begin Step 5 from
the updated parent branch.

## Why We Created This Integration Structure

PR #26 is not a small patch. It looks like a fork sync that brings in GPU
analysis, production-pipeline changes, DSSP refactoring, iPSAE tooling,
ColabFold conversion changes, metadata changes, Nextflow updates, and broad
documentation updates in one PR.

We decided not to merge it wholesale. Instead, we created a parent integration
branch and step branches so each area can be reviewed, verified, and merged
incrementally.

The governing plan is
[`docs/PR26_INTEGRATION_HANDOVER_PLAN.md`](./PR26_INTEGRATION_HANDOVER_PLAN.md).
That file remains the canonical checklist. This decision log explains the
context behind the completed steps and the current state.

## Branch And Commit Timeline

Initial setup:

- Fetched PR #26 to `review/pr-26`.
- Created parent integration branch `integration-pr-26-gpu` from the PR
  snapshot.
- Added the handover checklist:
  - `6e3c8dd docs: add PR 26 integration handover plan`

Step 1 branch:

- Branch: `integration-pr-26-gpu-step-1-dependencies`
- Commits:
  - `c94cd71 fix: make production dependencies optional at import`
  - `f96bbe6 fix: restore analysis metadata test target`
  - `ee769be docs: mark step 1 dependency work complete`
- Merged to parent:
  - `187c6d3 merge integration-pr-26-gpu-step-1-dependencies`

Step 2 branch:

- Branch: `integration-pr-26-gpu-step-2-hygiene`
- Commits reported by the Step 2 agent:
  - `48168ef chore: strip PR 26 whitespace noise`
  - `5eea333 docs: update step 2 hygiene status`
  - `717fddb docs: correct step 2 status commit note`
- Merged to parent:
  - `49483c6 merge integration-pr-26-gpu-step-2-hygiene`

Post Step 1/2 parent update:

- `ce1acb4 docs: update handover after merging steps 1 and 2`

Step 3 branch:

- Branch: `integration-pr-26-gpu-step-3-install-docs`
- Commits reported by the Step 3 agent:
  - `e2392e3 docs: align install instructions with dependency extras`
  - `ec10880 docs: update step 3 install docs status`
  - `add332b docs: clarify core uv install command`
- Merged to parent:
  - `41fa2f0 merge integration-pr-26-gpu-step-3-install-docs`

Step 4 branch:

- Branch: `integration-pr-26-gpu-step-4-dssp`
- Commits:
  - `4c6c347 fix: align DSSP algorithm defaults`
  - `6196ae1 test: cover internal DSSP algorithms`
  - `ffa6647 docs: update step 4 dssp status`
  - `8508d7a fix: restore mkdssp as DSSP default`
  - `a9a9d15 docs: clarify DSSP install requirement`
- Merge status:
  - Not merged into `integration-pr-26-gpu` at the time this decision log was
    written.

## Decisions Made

### 1. Integrate PR #26 Gradually

Decision:

- Treat PR #26 as an integration source branch, not a merge-ready PR.
- Use one step branch per integration area.
- Merge reviewed step branches back into `integration-pr-26-gpu` only after
  verification.

Reasoning:

- The PR is broad and fork-shaped.
- Smaller branches make regressions easier to isolate.
- Future chats can focus on one area at a time.

### 2. Preserve Core Install Without GPU/Production Dependencies

Decision:

- Core imports and CLI startup must not require production/GPU-only packages.
- `orjson` remains a core dependency because core code imports it.
- `biotite`, `pydssp`, `torch`, and `fastpdb` belong in the `production` extra.
- `torch_cluster` remains a separately documented install because its wheel must
  match the installed PyTorch/CUDA build.

Implemented in Step 1:

- Made DSSP and CIF/BCIF optional-package imports lazy where needed.
- Moved `biotite` and `pydssp` into the production extra.
- Removed unused `mdtraj`.
- Added/kept `orjson` in core requirements.
- Verified `main.py --help` works without production-only packages.

### 3. Restore Missing Test Target From Fork-Local Code

Decision:

- The copied test `tests/test_shard_analysis_metadata.py` must not import
  absent modules from `slurm-scaling/pipeline`.
- Bring the needed behavior into this repo under a normal package module.

Implemented in Step 1:

- Added `afdb_integration_kit.analysis_metadata`.
- Updated `tests/test_shard_analysis_metadata.py` to import that module.
- Full pytest passed after this change.

### 4. Mechanical Hygiene Is Separate From Behavior

Decision:

- PR-wide whitespace and end-of-file cleanup belongs in a mechanical Step 2
  branch.
- Do not reformat vendored code unless necessary.

Implemented in Step 2:

- Removed PR-introduced whitespace noise.
- Left vendored `afdb_integration_kit/ipsae/deps/json.hpp` untouched.

### 5. Documentation Must Match Dependency Split

Decision:

- Core install should be documented as `uv sync --locked --no-dev`.
- Contributors who need dev tools/tests can use `uv sync --locked`.
- Production install is `uv pip install '.[production]'`.
- `torch_cluster` is documented separately with a PyTorch/CUDA compatibility
  note.
- Docker remains core-only for Python dependencies unless a future step
  intentionally changes that.

Implemented in Step 3:

- Removed the missing `environment.yml` recommendation.
- Clarified core vs production vs GPU installs.
- Clarified Docker behavior.

### 6. Preserve Original DSSP Behavior By Default

Decision:

- The original DSSP behavior used external `mkdssp`.
- The shared library and CLI default must remain `mkdssp`.
- New Python algorithms are useful but must be opt-in:
  - `psea`
  - `pydssp`
  - `tmalign`
- The standalone production pipeline is allowed to default to `pydssp`, but it
  must be explicitly documented as production-specific.

Implemented in Step 4:

- Added `mkdssp` as a first-class algorithm option.
- Restored `DEFAULT_ALGORITHM = "mkdssp"`.
- Made `run_dssp(..., algorithm="mkdssp")` call the external subprocess path.
- Added CLI validation/help for all four algorithms.
- Added tests for default `mkdssp`, monkeypatched subprocess behavior, missing
  `mkdssp` handling, `tmalign` CIF output, and optional `psea`/`pydssp`.
- Updated README wording so DSSP install is required for the default
  `run-dssp`/`batch-dssp` path and for Nextflow workflows, while the production
  pipeline defaults to `pydssp`.

## Verification Results Reported So Far

After merging Steps 1 and 2 into the parent branch:

- `.venv/bin/pytest -q`: passed with `38 passed, 1 skipped`.
- `python3 -m compileall -q main.py afdb_integration_kit uniprot scripts tests`:
  passed.
- `git diff --check`: passed.
- `uv lock --locked`: passed with `UV_CACHE_DIR=/tmp/uv-cache`.
- `uv export --locked --no-hashes --output-file=requirements.txt --no-dev`:
  passed.
- `uv export --locked --no-hashes --extra production --output-file
  /tmp/production-req.txt --no-dev`: passed.

After Step 3:

- `uv sync --locked --no-dev --dry-run`: passed, with uv-reported local
  environment drift.
- `uv run main.py --help`: passed.
- `uv run main.py list-validations`: passed.
- `uv run python scripts/production_pipeline.py --help`: passed.
- `uv run python scripts/prepare_inputs.py --help`: passed.
- `.venv/bin/python -m pytest -q`: passed with `38 passed, 1 skipped`.
- `git diff --check`: passed.

After Step 4:

- `uv run pytest tests/test_dssp.py -q`: passed with `8 passed, 1 skipped`.
- `uv run pytest tests/test_pdb.py tests/validation/test_runner.py -q`: passed
  with `21 passed`.
- `uv run pytest -q`: reported by Step 4 agent as `46 passed, 2 skipped`.
- `.venv/bin/python -m compileall -q afdb_integration_kit/dssp tests`: passed.
- `git diff --check`: passed.

## Known Caveats And Follow-Up Needs

- Step 4 is complete but not yet merged into `integration-pr-26-gpu`.
- After merging Step 4, rerun the parent-branch verification before starting
  Step 5.
- The production pipeline intentionally defaults to `pydssp`; this is the one
  explicit DSSP-default exception. The shared CLI/library default is `mkdssp`.
- Python algorithms (`psea`, `pydssp`, `tmalign`) are 3-state approximations and
  should not be treated as byte-for-byte `mkdssp` replacements.
- `torch_cluster` is not in `pyproject.toml`; it is intentionally documented as
  a separate install due to PyTorch/CUDA wheel compatibility.
- Some future steps will require stronger review because they touch scientific
  correctness and output formats, especially ModelCIF, ColabFold conversion,
  GPU clash/interface analysis, iPSAE, and Nextflow.

## Immediate Next Action

Merge Step 4 back into the parent branch:

1. Switch to `integration-pr-26-gpu`.
2. Merge `integration-pr-26-gpu-step-4-dssp` with a normal merge commit.
3. Resolve conflicts if any, especially in:
   - `docs/PR26_INTEGRATION_HANDOVER_PLAN.md`
   - `README.md`
4. Run parent verification:
   - `git status --short --branch`
   - `git log --oneline --decorate --graph -14`
   - `git diff --check`
   - `.venv/bin/python -m pytest -q`
   - `.venv/bin/python -m compileall -q main.py afdb_integration_kit uniprot scripts tests`
   - `UV_CACHE_DIR=/tmp/uv-cache uv lock --locked`
   - `UV_CACHE_DIR=/tmp/uv-cache uv export --locked --no-hashes --output-file=requirements.txt --no-dev`
   - `UV_CACHE_DIR=/tmp/uv-cache uv export --locked --no-hashes --extra production --output-file /tmp/production-req.txt --no-dev`
5. Leave `integration-pr-26-gpu` clean.
6. Then Step 5 can start from the updated parent branch.

## Next Planned Step After Step 4 Merge

Step 5 in
[`docs/PR26_INTEGRATION_HANDOVER_PLAN.md`](./PR26_INTEGRATION_HANDOVER_PLAN.md)
is:

> CIF To BCIF Conversion Review

Recommended model for Step 5:

- Model: `GPT-5.4`
- Reasoning: `medium`

Purpose of Step 5:

- Review the new Biotite-first CIF to BCIF conversion path.
- Decide whether Biotite should remain optional/production-only or move back
  into core.
- Add tests for `.bcif`, `.bcif.gz`, fallback behavior, missing value masks, and
  temporary-file safety.
- Make sure the existing Mol* external path remains available and behavior is
  clear.

## Prompt Pattern For Future Coordinator Sessions

Future coordinator sessions should start by reading:

1. This decision log.
2. `docs/PR26_INTEGRATION_HANDOVER_PLAN.md`.

Then the coordinator should inspect branch state and return the exact next task
prompt plus model recommendation. This avoids stuffing every future worker
prompt with all historical context.
