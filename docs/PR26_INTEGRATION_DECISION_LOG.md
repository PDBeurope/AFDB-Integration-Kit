# PR 26 Integration Decision Log

This document summarizes the PR #26 integration work and decisions made so far.
It is intended as a compact handover for future Codex sessions. Use it together
with the detailed task checklist in
[`docs/PR26_INTEGRATION_HANDOVER_PLAN.md`](./PR26_INTEGRATION_HANDOVER_PLAN.md).

## Current State

- PR being integrated: GitHub PR #26 for `PDBeurope/AFDB-Integration-Kit`.
- PR snapshot branch: `review/pr-26`.
- Parent integration branch: `integration-pr-26-gpu`.
- Current verified parent branch:
  `integration-pr-26-gpu`.
- Step branches completed so far:
  - `integration-pr-26-gpu-step-1-dependencies`
  - `integration-pr-26-gpu-step-2-hygiene`
  - `integration-pr-26-gpu-step-3-install-docs`
  - `integration-pr-26-gpu-step-4-dssp`
  - `integration-pr-26-gpu-step-5-cif2bcif`
- Merged into `integration-pr-26-gpu` so far:
  - Step 1 dependencies/test baseline
  - Step 2 mechanical hygiene
  - Step 3 install documentation
  - Step 4 DSSP refactor branch
  - Step 5 CIF to BCIF conversion review

The immediate next action should be to create `integration-pr-26-gpu-step-6-colabfold-manifest`
from the verified parent branch and begin the Step 6 ColabFold/manifest
review there.

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
  - `3e97327 docs: add PR 26 integration decision log`
- Merge status:
  - Merged into `integration-pr-26-gpu` with merge commit `f2e5bb0`.

Step 5 branch:

- Branch: `integration-pr-26-gpu-step-5-cif2bcif`
- Commits:
  - `98c6600 fix: restore molstar as default cif2bcif backend`
  - `2778ca9 test: cover cif2bcif backends and temp safety`
  - `1488930 docs: record step 5 cif2bcif status`
- Merge status:
  - Merged into `integration-pr-26-gpu` with merge commit `2c6d225`.

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

### 7. Preserve Original CIF To BCIF Behavior By Default

Decision:

- The original toolkit converted CIF to BCIF by delegating to the external Mol*
  `cif2bcif` command.
- Step 5 should preserve that original behavior as the default/source-of-truth
  conversion path.
- The PR's Biotite in-process conversion logic should not be discarded, but it
  must be additive and explicit rather than silently replacing Mol*.
- Biotite should remain in the `production` extra and must continue to be lazily
  imported. Do not move Biotite into core dependencies for Step 5.
- Preferred backend shape:
  - `molstar`: use only external Mol* `cif2bcif`; this is the original and
    safest default behavior.
  - `biotite`: use only the Biotite in-process converter; useful for targeted
    testing and intentionally provisioned production environments.
  - `auto`: try Mol* first, then fall back to Biotite when Mol* is unavailable
    or fails.
- Do not make Biotite-first the implicit default.

Reasoning:

- Mol* `cif2bcif` is the established converter used by the original toolkit and
  is the lowest-risk path for downstream Mol*/gemmi/AFDB compatibility.
- The Biotite implementation may be useful in environments where Node/Mol* is
  hard to provision, but it makes this repository responsible for BinaryCIF
  encoding details such as column typing and `.`/`?` masks.
- Keeping Biotite explicit or fallback-only incorporates the PR's work without
  increasing the default integration risk.

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

After merging Step 4 into the parent branch:

- Merge commit: `f2e5bb0`.
- `uv run pytest tests/test_pdb.py tests/validation/test_runner.py -q`: passed
  with `21 passed`.
- `uv run pytest tests/test_dssp.py -q`: passed with `8 passed, 1 skipped, 1
  warning`.
- `.venv/bin/python -m compileall -q afdb_integration_kit/dssp tests`: passed.
- `git diff --check`: passed.
- `.venv/bin/python -m pytest -q`: passed with `46 passed, 2 skipped, 1
  warning`.

After Step 5 on `integration-pr-26-gpu-step-5-cif2bcif`:

- `uv run pytest tests/test_cif2bcif.py -q`: passed with `9 passed`.
- `uv run main.py run-cif2bcif --help`: passed.
- `uv run main.py batch-cif2bcif --help`: passed.
- `.venv/bin/python -m compileall -q afdb_integration_kit/cif2bcif tests`:
  passed.
- `git diff --check`: passed.
- `.venv/bin/python -m pytest -q`: passed with `55 passed, 2 skipped, 1
  warning`.

After merging Step 5 into the parent branch:

- Merge commit: `2c6d225`.
- `uv run pytest tests/test_cif2bcif.py -q`: passed with `9 passed`.
- `uv run main.py run-cif2bcif --help`: passed.
- `uv run main.py batch-cif2bcif --help`: passed.
- `.venv/bin/python -m compileall -q afdb_integration_kit/cif2bcif tests`:
  passed.
- `git diff --check`: passed.
- `.venv/bin/python -m pytest -q`: passed with `55 passed, 2 skipped, 1
  warning`.

## Known Caveats And Follow-Up Needs

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

Continue on `integration-pr-26-gpu-step-6-colabfold-manifest`, which was
created from the verified `integration-pr-26-gpu` parent branch, and begin the
Step 6 review from
[`docs/PR26_INTEGRATION_HANDOVER_PLAN.md`](./PR26_INTEGRATION_HANDOVER_PLAN.md).

## Next Planned Step After Step 5

Step 6 in
[`docs/PR26_INTEGRATION_HANDOVER_PLAN.md`](./PR26_INTEGRATION_HANDOVER_PLAN.md)
is:

> ColabFold Converter And Manifest Resolver Review

Recommended model for Step 6:

- Model: `GPT-5.4`
- Reasoning: `medium`

Purpose of Step 6:

- Review `afdb_integration_kit/colabfold/converter.py` and
  `afdb_integration_kit/manifest/resolver.py`.
- Confirm chain-span parsing, PAE rounding, caching behavior, and AFID
  normalization assumptions before merging.

## Prompt Pattern For Future Coordinator Sessions

Future coordinator sessions should start by reading:

1. This decision log.
2. `docs/PR26_INTEGRATION_HANDOVER_PLAN.md`.

Then the coordinator should inspect branch state and return the exact next task
prompt plus model recommendation. This avoids stuffing every future worker
prompt with all historical context.
