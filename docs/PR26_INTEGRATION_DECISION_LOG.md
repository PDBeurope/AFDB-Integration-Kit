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
  - `integration-pr-26-gpu-step-6-colabfold-manifest`
  - `integration-pr-26-gpu-step-7-gpu-analysis`
- Merged into `integration-pr-26-gpu` so far:
  - Step 1 dependencies/test baseline
  - Step 2 mechanical hygiene
  - Step 3 install documentation
  - Step 4 DSSP refactor branch
  - Step 5 CIF to BCIF conversion review
  - Step 6 ColabFold converter and manifest resolver review
  - Step 7 GPU clash/interface analysis package review

The current in-progress integration branch is
`integration-pr-26-gpu-step-8-ipsae`. Step 8 has been implemented and verified
on that branch, but it has not been merged back into `integration-pr-26-gpu`
yet.

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

Step 6 branch:

- Branch: `integration-pr-26-gpu-step-6-colabfold-manifest`
- Commits:
  - `830e12f docs: mark step 6 branch prepared`
  - `a7d28c6 fix: preserve colabfold manifest semantics`
  - `de7e51f docs: record step 6 colabfold review status`
  - `77ccadc test: add curated colabfold real example fixtures`
  - `c7549eb test: normalize colabfold real fixtures`
  - `30de578 docs: update step 6 commit list`
  - `0024936 docs: document colabfold fixture handoff`
- Merge status:
  - Merged into `integration-pr-26-gpu` with merge commit `4c4a158`.

Step 7 branch:

- Branch: `integration-pr-26-gpu-step-7-gpu-analysis`
- Scope note before implementation:
  - `e9fa837 docs: clarify step 7 gpu analysis scope`
- Current branch status:
  - Package import is now lazy: `import afdb_integration_kit.gpu` no longer
    requires `torch`, `fastpdb`, or `biotite`.
  - Production analysis entry points now support `device="auto"` and fail
    early with a clear error when `device="cuda"` is requested without CUDA.
  - `parse.py` now imports `fastpdb` and Biotite lazily, so parsing
    dependencies are only required when PDB parsing is actually invoked.
  - Added `tests/test_gpu_analysis.py` for import behavior, dependency error
    messaging, schema conversion, device resolution, and CPU/fallback coverage
    when PyTorch is available.
- Merge status:
  - Merged into `integration-pr-26-gpu` with merge commit `0529ffa`.

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

### 8. Preserve Global Chain Spans In ColabFold JSON Outputs

Decision:

- Keep `chains.sequenceStart` and `chains.sequenceEnd` aligned to the global
  flattened pLDDT/PAE indexing used by the original toolkit outputs.
- Keep the gemmi parser as an implementation detail, but do not change the
  external chain-span contract.

Reasoning:

- `origin/main` emitted global chain spans.
- The PR branch changed those spans to per-chain local numbering while still
  emitting a global `residueNumber` array, which creates a mixed indexing
  contract.
- Step 6 keeps the performance improvement and parser fallback, but restores
  the original output semantics.

### 9. Do Not Guess Ambiguous UniProt Accessions

Decision:

- If an AF ID maps to multiple candidate accessions and no DuckDB-plus-meta
  evidence is available to disambiguate them, fail that AF ID instead of
  choosing the first accession alphabetically.
- Keep the pLDDT-length-based disambiguation path when the necessary evidence
  is available.

Reasoning:

- Alphabetical fallback is not scientifically meaningful.
- A failed mapping is easier to detect and repair than a silently incorrect
  accession assignment.
- This keeps Step 6 conservative while still allowing deterministic resolution
  when the accession lengths and ColabFold metadata provide direct evidence.

### 10. Use Single AF-Style IDs For Real ColabFold Fixtures

Decision:

- Curated real ColabFold fixtures should use single AF-style `model_entity_id`
  names in directories, copied files, and fixture metadata.
- Heterodimer source files that were named with two component AF IDs should be
  normalized to a single fixture AF ID; the original composite source ID and
  component reassigned AF IDs are retained in `manifest.json`.
- Chain UniProt accessions should be recorded in `manifest.json` so later
  integration and end-to-end tests can build real-style chain manifests without
  relying on filename parsing.

Reasoning:

- The production `model_entity_id` convention is a single AF-style ID, not a
  composite ID made from both components.
- Keeping source IDs as metadata preserves provenance while making tests target
  the naming contract the toolkit is expected to support.

### 11. Keep The GPU Package Optional At Import Time

Decision:

- `afdb_integration_kit.gpu` must remain importable in a core-only
  environment.
- The package `__init__` should export production analysis APIs lazily instead
  of importing Torch, parsing, and GPU modules eagerly.
- Production analysis modules should raise clear, actionable
  optional-dependency errors when `torch`, `fastpdb`, or related packages are
  missing.
- `parse.py` should load `fastpdb` and Biotite lazily so
  `analyze_proteins()`/schema helpers are not coupled to PDB parsing.

Implemented in Step 7:

- Added `afdb_integration_kit.gpu._runtime` helpers for optional dependency
  imports and device resolution.
- Reworked `afdb_integration_kit.gpu.__init__` to lazy-load production
  symbols.
- Verified `import afdb_integration_kit.gpu`, `import
  afdb_integration_kit.gpu.parse`, and `main.py --help` work in this sandbox
  even though `.venv` does not contain `torch` or `fastpdb`.

### 12. Make Device Selection Explicit And Conservative

Decision:

- Public GPU analysis entry points should accept `device="cpu"`,
  `device="cuda"`, and `device="auto"`.
- `device="auto"` should resolve to CUDA when available, otherwise CPU.
- `device="cuda"` should fail before parsing or batching work starts if CUDA is
  unavailable.
- CPU execution is supported for small correctness workloads; CUDA remains the
  intended production path.

Implemented in Step 7:

- Updated the public analysis and batching entry points to normalize device
  selection through one helper.
- Changed the CLI default in `afdb_integration_kit.gpu.analyze` to
  `device="auto"`.
- Documented the CPU-only verification boundary: API/correctness checks were
  exercised here, but CUDA throughput and GPU memory behavior were not.

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

After Step 6 on `integration-pr-26-gpu-step-6-colabfold-manifest`:

- `UV_CACHE_DIR=/tmp/uv-cache uv run pytest tests/test_colabfold_converter.py
  -q`: passed with `14 passed`.
- `UV_CACHE_DIR=/tmp/uv-cache uv run pytest tests/test_manifest_resolver.py
  -q`: passed with `4 passed`.
- `.venv/bin/python -m compileall -q afdb_integration_kit/colabfold
  afdb_integration_kit/manifest tests`: passed.
- `git diff --check`: passed.
- `.venv/bin/python -m pytest -q`: passed with `73 passed, 1 skipped, 1
  warning`.
- Real-fixture coverage now includes all 9 curated ColabFold examples under
  `tests/fixtures/colabfold_real_examples`.

After merging Step 6 into the parent branch:

- Merge commit: `4c4a158`.
- `UV_CACHE_DIR=/tmp/uv-cache uv run pytest tests/test_colabfold_converter.py
  -q`: passed with `14 passed`.
- `UV_CACHE_DIR=/tmp/uv-cache uv run pytest tests/test_manifest_resolver.py
  -q`: passed with `4 passed`.
- `.venv/bin/python -m compileall -q afdb_integration_kit/colabfold
  afdb_integration_kit/manifest tests`: passed.
- `git diff --check`: passed.
- `.venv/bin/python -m pytest -q`: passed with `73 passed, 1 skipped, 1
  warning`.

After Step 7 on `integration-pr-26-gpu-step-7-gpu-analysis`:

- `.venv/bin/python -m pytest -q tests/test_gpu_analysis.py`: passed with
  `6 passed, 1 skipped`.
- `.venv/bin/python -m pytest -q tests/test_cif2bcif.py tests/test_dssp.py
  tests/test_gpu_analysis.py`: passed with `23 passed, 2 skipped, 1 warning`.
- `.venv/bin/python -m compileall -q afdb_integration_kit/gpu
  tests/test_gpu_analysis.py`: passed.
- `.venv/bin/python main.py --help`: passed.
- `.venv/bin/python scripts/production_pipeline.py --help`: passed.
- Import probe in `.venv`:
  `afdb_integration_kit.gpu`, `.gpu.protein`, `.gpu.parse`, and `.gpu.schema`
  all import successfully without production dependencies installed.
- Environment note: `.venv` is missing `torch`, `fastpdb`, and
  `torch_cluster`, so the Step 7 CPU execution/fallback test is skip-marked and
  CUDA behavior remains unverified here.

After merging Step 7 into the parent branch:

- Merge commit: `0529ffa`.
- `git diff --check`: passed.
- `.venv/bin/python -m compileall -q afdb_integration_kit/gpu tests`: passed.
- `.venv/bin/python -m pytest -q tests/test_gpu_analysis.py`: passed with
  `6 passed, 1 skipped`.
- `.venv/bin/python -m pytest -q tests/test_cif2bcif.py tests/test_dssp.py
  tests/test_gpu_analysis.py`: passed with `23 passed, 2 skipped, 1 warning`.
- `.venv/bin/python main.py --help`: passed.
- `.venv/bin/python scripts/production_pipeline.py --help`: passed.
- `.venv/bin/python -m pytest -q`: passed with `79 passed, 2 skipped, 1
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

Review Step 8 on `integration-pr-26-gpu-step-8-ipsae` for merge readiness, or
continue with the next scoped branch only after Step 8 is accepted. Do not
start Step 9 work on this branch.

## Step 8 Plan Reference

Step 8 in
[`docs/PR26_INTEGRATION_HANDOVER_PLAN.md`](./PR26_INTEGRATION_HANDOVER_PLAN.md)
is:

> iPSAE C++ Tool Review

Recommended model for Step 8:

- Model: `GPT-5.4`
- Reasoning: `medium`

Purpose of Step 8:

- Review `afdb_integration_kit/ipsae/ipsae_cpp.cpp` and its build path.
- Decide whether vendored `json.hpp` is acceptable and document
  source/version/license if retained.
- Confirm the C++ tool builds locally where compiler support is available.
- Add a minimal fixture and smoke test for expected CSV output where feasible.

## Step 8 Outcome

Branch:

- `integration-pr-26-gpu-step-8-ipsae`

Scope completed on the branch:

- Reviewed `afdb_integration_kit/ipsae/ipsae_cpp.cpp`,
  `afdb_integration_kit/ipsae/Makefile`,
  `afdb_integration_kit/ipsae/deps/json.hpp`,
  `afdb_integration_kit/ipsae/README.md`,
  `scripts/production_pipeline.py`, and the repo/script references that mention
  iPSAE.
- Kept `deps/json.hpp` vendored and documented its provenance/license in the
  local iPSAE README.
- Updated the iPSAE `Makefile` so the default `make` path prefers an already
  installed Eigen tree before attempting a network fetch into `deps/`.
- Added `tests/test_ipsae_cpp.py` as a minimal build-and-run smoke test for the
  batch CSV contract.

Step 8 decisions:

- Vendored `nlohmann/json` is acceptable here.
  - Version in `deps/json.hpp`: `3.11.3`.
  - License in the retained upstream header: `MIT`.
  - Reasoning: it is header-only, already carries SPDX metadata, and keeping it
    in-repo avoids introducing a separate package-manager dependency for the
    production pipeline's C++ build path.
- Eigen remains non-vendored.
  - Reasoning: the repo already had a Makefile flow that fetched Eigen on
    demand; Step 8 kept that model but made it less network-dependent by
    preferring system/local Eigen headers when present.
- The Step 8 code change stayed scoped to build correctness and testability.
  - No scoring algorithm rewrite was attempted.

Numerical/scientific caveats recorded during Step 8 review:

- `pae_cutoff` directly controls which residue pairs contribute to ipSAE, but
  `dist_cutoff` is only used for the reported `dist_nres1`/`dist_nres2` counts.
  It does not alter the ipSAE score itself.
- pDockQ uses a fixed 8.0 A contact cutoff in the C++ implementation.
- LIS uses a fixed 12.0 A PAE cutoff in the C++ implementation.
- Missing `CB` atoms fall back to `CA` coordinates; the code comment now states
  this explicitly rather than implying it is glycine-only.
- The batch summary column naming assumes single-character chain IDs (for
  example `ipsae_AB` and `ipsae_BA`).
- `uniprot/scripts/batch_ipsae.py` appears fork-stale relative to the reviewed
  C++ binary path and CLI. Step 8 left it unchanged because the active
  production path uses `scripts/production_pipeline.py`; if that UniProt script
  is still meant to be supported, fix or retire it during Step 10.

Step 8 verification on the branch:

- `make check` in `afdb_integration_kit/ipsae`: passed, resolving Eigen from
  `/usr/include/eigen3`.
- `make clean && make` in `afdb_integration_kit/ipsae`: passed and built
  `ipsae_cpp`.
  - Caveat: the static link emitted a `libgomp.a`/`dlopen` warning from the
    host toolchain, but the build completed successfully.
- `.venv/bin/python -m pytest -q tests/test_ipsae_cpp.py
  tests/test_shard_analysis_metadata.py`: passed with `3 passed`.
- `.venv/bin/python -m compileall -q afdb_integration_kit tests`: passed.
- `git diff --check`: passed.
- `.venv/bin/python -m pytest -q`: passed with `80 passed, 2 skipped, 1 warning`.

## Prompt Pattern For Future Coordinator Sessions

Future coordinator sessions should start by reading:

1. This decision log.
2. `docs/PR26_INTEGRATION_HANDOVER_PLAN.md`.

Then the coordinator should inspect branch state and return the exact next task
prompt plus model recommendation. This avoids stuffing every future worker
prompt with all historical context.
