# PR 26 GPU Integration Handover Plan

This plan is for gradually integrating GitHub PR #26 into the AFDB Integration
Kit without merging the fork wholesale. Future chats should continue from the
branch and step named below, run the listed verification commands, and update
this checklist as work moves from pending to complete.

Current branch structure:

- Parent integration branch: `integration-pr-26-gpu`
- PR snapshot branch: `review/pr-26`
- Base branch for comparison: `origin/main`
- Step branches should be named `integration-pr-26-gpu-step-N-<short-topic>`

General rules for every step:

- [ ] (Model: GPT-5.4) Start by checking `git status --short --branch` and
  `git log --oneline --decorate -5`.
- [ ] (Model: GPT-5.4) Preserve user or prior-agent changes; do not revert
  unrelated edits.
- [ ] (Model: GPT-5.4) Keep each step on its own sub-branch from the parent
  integration branch unless the user explicitly asks to merge it back.
- [ ] (Model: GPT-5.4) Prefer small commits grouped by behavior: packaging,
  code behavior, tests, docs, and mechanical cleanup should not be mixed unless
  the diff is trivial.
- [ ] (Model: GPT-5.4-mini) Run `git diff --check` before each final response.
- [ ] (Model: GPT-5.4-mini) Run the narrowest relevant tests first, then
  `pytest -q` before marking a step complete.
- [ ] (Model: GPT-5.4) If `uv lock` or `uv sync` needs network and the current
  sandbox has restricted network, use `uv lock --locked`, `uv export --locked`,
  and dry-run checks. State clearly when a full online lock refresh was not run.
- [ ] (Model: GPT-5.5) Use a code-review pass for each step before merging any
  step branch back into `integration-pr-26-gpu`.

## Step 1: Create Internal Branch And Fix Packaging/Test Baseline

- [x] (Model: GPT-5.4) Fetch PR #26 locally.
  - Evidence: `review/pr-26` points at commit `4cdcfc0`.
- [x] (Model: GPT-5.4) Create the parent integration branch from PR #26.
  - Branch created: `integration-pr-26-gpu`.
  - Note: `integration/pr-26-gpu` could not be created in this sandbox because
    Git could not create the nested ref path, so the flat branch name was used.
- [x] (Model: GPT-5.4) Review the initial failure state.
  - `pytest -q` originally failed during collection because
    `tests/test_shard_analysis_metadata.py` imported absent fork-local modules
    from `slurm-scaling/pipeline`.
  - Importing `main` originally required production-only packages such as
    `pydssp`.
  - `requirements.txt` did not include all new core runtime packages.
- [x] (Model: GPT-5.4) Move production/GPU-only imports behind runtime
  boundaries.
  - `afdb_integration_kit.dssp.dssp` no longer imports `pydssp` or `torch` at
    module import time.
  - `afdb_integration_kit.cif2bcif.convert` no longer imports Biotite at module
    import time.
  - Expected behavior: `main.py --help` works without production/GPU packages.
- [x] (Model: GPT-5.4) Adjust dependency metadata.
  - Keep `orjson` as a core dependency because multiple core modules import it.
  - Move `biotite` and `pydssp` into the `production` extra.
  - Remove unused `mdtraj` from project dependencies.
  - Keep `torch` and `fastpdb` in the `production` extra.
  - Leave `torch_cluster` as a separately documented install because it is not
    in `pyproject.toml`.
- [x] (Model: GPT-5.4) Replace missing external analysis-metadata imports with
  a repo-owned module.
  - Add `afdb_integration_kit.analysis_metadata`.
  - Update `tests/test_shard_analysis_metadata.py` to import that module
    directly.
- [x] (Model: GPT-5.4-mini) Verify Step 1.
  - `pytest -q`: expected `38 passed, 1 skipped`.
  - `python -m compileall -q main.py afdb_integration_kit uniprot scripts tests`:
    expected pass.
  - targeted `flake8` on touched files: expected pass.
  - `git diff --check`: expected pass.
  - `uv lock --locked`: expected pass.
  - `uv export --locked --no-hashes --output-file=requirements.txt --no-dev`:
    expected pass.
  - `uv export --locked --no-hashes --extra production --output-file
    /tmp/production-req.txt --no-dev`: expected pass and should include
    `biotite`, `pydssp`, `torch`, and `fastpdb`, but not `mdtraj`.
  - `uv sync --locked --no-dev --dry-run`: expected to resolve. In this
    sandbox it reports local environment drift, including replacing
    `orjson==3.11.5` with locked `orjson==3.11.4`.
- [x] (Model: GPT-5.4) Commit only this handover plan on the parent integration
  branch before committing Step 1 implementation changes.
- [x] (Model: GPT-5.4) Create Step 1 sub-branch:
  `integration-pr-26-gpu-step-1-dependencies`.
- [x] (Model: GPT-5.4) Commit the Step 1 implementation changes on the Step 1
  sub-branch.
  - Commit 1: `c94cd71 fix: make production dependencies optional at import`
  - Commit 2: `f96bbe6 fix: restore analysis metadata test target`
  - Commit 3: this checklist/status update.

## Step 2: Mechanical Hygiene And Repository Fit

- [x] (Model: GPT-5.4-mini) Create branch
  `integration-pr-26-gpu-step-2-hygiene` from `integration-pr-26-gpu` after
  Step 1 is reviewed or merged back.
- [x] (Model: GPT-5.4-mini) Run `git diff --check origin/main...HEAD` and list
  all whitespace issues introduced by PR #26.
  - Before cleanup, the issues were concentrated in:
    `afdb_integration_kit/colabfold/converter.py`,
    `afdb_integration_kit/gpu/README.md`,
    `afdb_integration_kit/ipsae/ipsae_cpp.cpp`,
    `afdb_integration_kit/validation/validators/_parallel.py`,
    `docs/PIPELINE_OPTIMIZATION_ANALYSIS.md`, `main.py`,
    `scripts/production_pipeline.py`, `tests/fixtures/setup_mock_data.py`,
    `uniprot/scripts/batch_export_modelcif_input.py`,
    `uniprot/scripts/batch_ipsae.py`, `uniprot/scripts/batch_validate_assets.py`,
    and `workflow/end_to_end_with_validation_multibatch.nf`.
- [x] (Model: GPT-5.4-mini) Fix trailing whitespace and end-of-file issues
  across the PR.
  - Do not reformat vendored `afdb_integration_kit/ipsae/deps/json.hpp` unless
    necessary.
  - Prefer mechanical whitespace-only commits.
- [x] (Model: GPT-5.4) Run pre-commit-equivalent checks available locally:
  - `git diff --check`
  - `.venv/bin/flake8 <changed-python-files> --max-line-length=88`
  - `python -m compileall -q main.py afdb_integration_kit uniprot scripts tests`
  - `flake8` still reports existing style issues in the touched Python files;
    this step stayed scoped to whitespace-only cleanup.
  - `python3 -m compileall -q main.py afdb_integration_kit uniprot scripts tests`
    passed in this sandbox because `python` is not on `PATH`.
- [x] (Model: GPT-5.4) Decide whether to add a repo-level formatter config or
  limit the change to cleaning PR-introduced whitespace.
  - Decision: limit to PR-introduced whitespace cleanup; no formatter config
    added.
- [x] (Model: GPT-5.4) Commit mechanical cleanup separately from behavior.
  - Cleanup commit: `48168ef chore: strip PR 26 whitespace noise`
  - Step 2 plan/status update commit: `5eea333 docs: update step 2 hygiene status`

Note: Step 1 and Step 2 have been merged back into `integration-pr-26-gpu`.
After the merge, `pytest -q` is expected to use the repo-owned
`afdb_integration_kit.analysis_metadata` module and the optional production
dependency import boundaries from Step 1.

## Step 3: Documentation And Install Instructions

- [x] (Model: GPT-5.4) Create branch
  `integration-pr-26-gpu-step-3-install-docs`.
- [x] (Model: GPT-5.4) Audit README install sections for claims introduced by
  PR #26.
  - Decision: pip/uv is the intended path for this integration step.
  - Removed the missing `environment.yml` recommendation instead of adding an
    untested conda environment.
- [x] (Model: GPT-5.4) Align install docs with dependency metadata.
  - Core install should be enough for CLI help, validation, metadata, UniProt,
    ColabFold conversion, ModelCIF/PDB generation, and non-production scripts.
  - Production install should mention `uv pip install '.[production]'`.
  - GPU install should mention `torch_cluster` separately with the correct
    PyTorch/CUDA compatibility note.
- [x] (Model: GPT-5.4) Check Dockerfile expectations.
  - Decision: Docker remains core-only for Python dependencies.
  - Kept the Dockerfile `requirements.txt` install unchanged and documented
    that the image does not install `.[production]` or `torch_cluster`.
- [x] (Model: GPT-5.4-mini) Verify all documented commands that can run without
  external datasets.
  - `uv sync --locked --no-dev --dry-run`: passed, with local environment drift
    noted by uv.
  - `uv run main.py --help`: passed.
  - `uv run main.py list-validations`: passed.
  - `uv run python scripts/production_pipeline.py --help`: passed.
  - `uv run python scripts/prepare_inputs.py --help`: passed.
  - `uv run main.py test`: fails in this sandbox because `cif2bcif` is not on
    `PATH`; README now documents it as an optional external toolchain check.
  - `.venv/bin/python -m pytest -q`: expected `38 passed, 1 skipped`.
  - `git diff --check`: expected pass.
- [x] (Model: GPT-5.4) Commit docs separately from Docker or build-system
  changes.
  - No Dockerfile or build-system changes were needed.

## Step 4: DSSP Refactor Review

- [x] (Model: GPT-5.4) Create branch
  `integration-pr-26-gpu-step-4-dssp`.
- [x] (Model: GPT-5.4) Review the replacement of external `mkdssp` subprocess
  behavior with internal algorithms in `afdb_integration_kit/dssp/dssp.py`.
- [x] (Model: GPT-5.4) Confirm desired default algorithm.
  - README currently indicates `pydssp` in production-pipeline contexts.
  - CLI `batch-dssp` default may differ; reconcile defaults intentionally.
  - Follow-up correction: `mkdssp` remains the shared default to preserve the
    historical CLI and library behavior. `pydssp` remains explicit only in the
    production pipeline, where Python production dependencies are documented.
  - `tmalign` remains the core-dependency algorithm for environments without
    production extras and is covered directly by tests.
- [x] (Model: GPT-5.4) Add tests for at least `tmalign`, because it has no
  production-only dependencies.
- [x] (Model: GPT-5.4) Add skip-marked tests for `psea` and `pydssp` when
  optional packages are not installed.
- [x] (Model: GPT-5.4) Validate output CIF structure.
  - Confirm `_struct_conf` and `_struct_conf_type` are generated correctly.
  - Confirm existing CIF categories are preserved.
  - Confirm files with no polymer residues fail cleanly.
- [x] (Model: GPT-5.4-mini) Run:
  - `pytest tests/test_pdb.py tests/validation/test_runner.py -q`
  - any new DSSP tests
  - `python -m compileall -q afdb_integration_kit/dssp tests`
  - Evidence: `uv run pytest tests/test_pdb.py tests/validation/test_runner.py
    -q` passed with `21 passed`.
  - Evidence: `uv run pytest tests/test_dssp.py -q` passed with `8 passed, 1
    skipped`; the skip was the real `mkdssp` executable check in this sandbox.
  - Evidence: `python -m compileall -q afdb_integration_kit/dssp tests` could
    not run because `python` is not on `PATH`; `.venv/bin/python -m compileall
    -q afdb_integration_kit/dssp tests` passed.
  - Evidence: `git diff --check` passed.
- [x] (Model: GPT-5.5) Review biological correctness risks before merging.
  - Risk note: `pydssp`, `psea`, and `tmalign` are all 3-state approximations
    rather than byte-for-byte `mkdssp` replacements. The shared default remains
    `mkdssp`; the Python algorithms should be treated as opt-in alternatives.

Note: Step 4 has been merged back into `integration-pr-26-gpu` with merge
commit `f2e5bb0`. Parent-branch verification after the merge:

- `uv run pytest tests/test_pdb.py tests/validation/test_runner.py -q`: passed
  with `21 passed`.
- `uv run pytest tests/test_dssp.py -q`: passed with `8 passed, 1 skipped, 1
  warning`.
- `.venv/bin/python -m compileall -q afdb_integration_kit/dssp tests`: passed.
- `git diff --check`: passed.
- `.venv/bin/python -m pytest -q`: passed with `46 passed, 2 skipped, 1
  warning`.

## Step 5: CIF To BCIF Conversion Review

- [x] (Model: GPT-5.4) Create branch
  `integration-pr-26-gpu-step-5-cif2bcif`.
- [ ] (Model: GPT-5.4) Review the new Biotite conversion path in
  `afdb_integration_kit/cif2bcif/convert.py`.
- [x] (Model: GPT-5.4) Decide the integration direction for Biotite.
  - Decision: preserve the original toolkit contract that Mol* `cif2bcif` is
    the default/source-of-truth conversion path.
  - Decision: do not move Biotite into core dependencies. Keep it in the
    `production` extra and import it lazily.
  - Decision: keep the PR's Biotite conversion work, but make it explicit and
    additive rather than silently replacing the original behavior.
- [x] (Model: GPT-5.4) Implement backend selection conservatively.
  - Add or preserve a backend option with values equivalent to:
    `molstar`, `biotite`, and `auto`.
  - `molstar`: run only the external Mol* `cif2bcif` command. This is the
    original behavior and should remain the safest/default behavior unless the
    user explicitly chooses otherwise.
  - `biotite`: run only the Biotite in-process converter. This is useful for
    targeted testing or production environments that intentionally install
    Biotite.
  - `auto`: try Mol* first, then fall back to Biotite if Mol* is unavailable or
    fails.
  - Do not make Biotite-first the implicit default.
- [x] (Model: GPT-5.4) Add tests for:
  - default Mol* command selection and subprocess behavior.
  - explicit Biotite backend behavior, skip-marked if Biotite is not installed.
  - `auto` fallback from Mol* to Biotite when Mol* is unavailable or fails.
  - no import-time Biotite requirement.
  - `.bcif` output creation.
  - `.bcif.gz` output creation.
  - missing value mask behavior for `.` and `?` in the Biotite backend.
- [x] (Model: GPT-5.4) Check temporary-file handling.
  - Ensure concurrent workers cannot collide on the same temp filename.
  - Ensure cross-device rename fallback is safe.
- [x] (Model: GPT-5.4) Update README/help text only as needed so users know
  Mol* is the default path and Biotite is optional/explicit/fallback.
- [x] (Model: GPT-5.4-mini) Run existing and new CIF/BCIF tests.
  - Evidence: `uv run pytest tests/test_cif2bcif.py -q` passed with `9 passed`.
  - Evidence: `uv run main.py run-cif2bcif --help` passed and documents
    `molstar` as the default backend with `biotite` and `auto` available.
  - Evidence: `uv run main.py batch-cif2bcif --help` passed and documents the
    same backend choices.
  - Evidence: `.venv/bin/python -m compileall -q afdb_integration_kit/cif2bcif
    tests` passed.
  - Evidence: `git diff --check` passed.
  - Evidence: `.venv/bin/python -m pytest -q` passed with `55 passed, 2
    skipped, 1 warning`.
- [x] (Model: GPT-5.5) Review compatibility with downstream Mol* and gemmi
  consumers before merging.
  - Risk note: keeping Mol* as the default backend preserves the original
    converter contract for downstream Mol*/gemmi consumers. The Biotite backend
    remains explicit or fallback-only so BinaryCIF encoding changes do not
    silently replace the established path.

Note: Step 5 has been merged back into `integration-pr-26-gpu` with merge
commit `2c6d225`. Parent-branch verification after the merge:

- `uv run pytest tests/test_cif2bcif.py -q`: passed with `9 passed`.
- `uv run main.py run-cif2bcif --help`: passed.
- `uv run main.py batch-cif2bcif --help`: passed.
- `.venv/bin/python -m compileall -q afdb_integration_kit/cif2bcif tests`:
  passed.
- `git diff --check`: passed.
- `.venv/bin/python -m pytest -q`: passed with `55 passed, 2 skipped, 1
  warning`.

## Step 6: ColabFold Converter And Manifest Resolver Review

- [ ] (Model: GPT-5.4) Create branch
  `integration-pr-26-gpu-step-6-colabfold-manifest`.
- [ ] (Model: GPT-5.4) Review `afdb_integration_kit/colabfold/converter.py`
  changes.
  - Confirm gemmi parsing gives identical chain spans to the previous parser
    for representative PDBs.
  - Confirm PAE rounding and `orjson.OPT_SERIALIZE_NUMPY` are intended.
  - Confirm DuckDB prefetch cache is safe across repeated conversions.
- [ ] (Model: GPT-5.4) Review `afdb_integration_kit/manifest/resolver.py`.
  - Confirm supported model ID formats: `AF-...`, `AF_...`, and
    `AF_..._AF_...`.
  - Confirm 16-digit AF IDs are accepted if that is the production convention.
  - Confirm hyphen/underscore normalization is consistent with filenames,
    manifests, and metadata.
  - Confirm ambiguous accession deduplication is scientifically valid.
- [ ] (Model: GPT-5.4) Add tests for model-ID classification and manifest row
  building.
- [ ] (Model: GPT-5.4) Add tests for chain mapping behavior in
  `convert_file`.
- [ ] (Model: GPT-5.4-mini) Run:
  - `pytest tests/test_colabfold_converter.py -q`
  - new manifest resolver tests
- [ ] (Model: GPT-5.5) Review data-model assumptions before merging.

## Step 7: GPU Clash/Interface Package Review

- [ ] (Model: GPT-5.4) Create branch
  `integration-pr-26-gpu-step-7-gpu-analysis`.
- [ ] (Model: GPT-5.4) Review `afdb_integration_kit/gpu/*` as a package.
  - Confirm imports are optional and do not affect core package import.
  - Consider making `afdb_integration_kit.gpu.__init__` lazy or minimal so
    missing `fastpdb`, `torch`, or `torch_cluster` produce clear errors.
- [ ] (Model: GPT-5.4) Check licensing headers and provenance for copied GPU
  code.
- [ ] (Model: GPT-5.4) Add CPU-compatible tests where possible.
  - Test schema conversion.
  - Test small clash/interface calculations if torch is available.
  - Skip GPU-specific tests unless CUDA is available.
- [ ] (Model: GPT-5.4) Check `torch_cluster` fallback behavior in
  `clashes.py`.
- [ ] (Model: GPT-5.5) Review performance-sensitive code for memory pressure
  and batching assumptions.

## Step 8: iPSAE C++ Tool Review

- [ ] (Model: GPT-5.4) Create branch
  `integration-pr-26-gpu-step-8-ipsae`.
- [ ] (Model: GPT-5.4) Review `afdb_integration_kit/ipsae/ipsae_cpp.cpp` and
  `Makefile`.
- [ ] (Model: GPT-5.4) Decide whether vendored `json.hpp` is acceptable.
  - If vendored, document source/version/license.
  - If not, switch to a system/package dependency.
- [ ] (Model: GPT-5.4) Confirm the C++ tool builds locally with available
  compilers.
- [ ] (Model: GPT-5.4) Add a minimal fixture and test for expected CSV output,
  if build tooling is available.
- [ ] (Model: GPT-5.4-mini) Run C++ build and smoke test where supported.
- [ ] (Model: GPT-5.5) Review numerical correctness and threshold handling.

## Step 9: Production Pipeline Scripts Review

- [ ] (Model: GPT-5.4) Create branch
  `integration-pr-26-gpu-step-9-production-pipeline`.
- [ ] (Model: GPT-5.4) Review `scripts/prepare_inputs.py` and
  `scripts/production_pipeline.py`.
- [ ] (Model: GPT-5.4) Confirm all paths and generated config files line up
  with the repo-owned modules.
- [ ] (Model: GPT-5.4) Confirm production mode does not require network/API.
- [ ] (Model: GPT-5.4) Confirm dev mode handles API failures and missing
  accessions predictably.
- [ ] (Model: GPT-5.4) Add tests for:
  - matched PDB/JSON pair detection.
  - symlink/copy behavior.
  - config generation.
  - resume/caching behavior if feasible with small fixtures.
- [ ] (Model: GPT-5.4) Decide whether `analysis_metadata.py` should remain a
  library module or move under `scripts/`.
- [ ] (Model: GPT-5.5) Review operational safety before merging: idempotency,
  cleanup, retry behavior, and failure reports.

## Step 10: UniProt Script And Template Review

- [ ] (Model: GPT-5.4) Create branch
  `integration-pr-26-gpu-step-10-uniprot`.
- [ ] (Model: GPT-5.4) Review changes under `uniprot/scripts/` and
  `afdb_integration_kit/uniprot/api.py`.
- [ ] (Model: GPT-5.4) Confirm REST API fetching is only used in dev/small-scale
  mode and not production mode.
- [ ] (Model: GPT-5.4) Confirm generated DuckDB schema matches existing
  extraction tooling.
- [ ] (Model: GPT-5.4) Restore or replace deleted UniProt tests if they still
  cover relevant behavior.
- [ ] (Model: GPT-5.4) Add tests for parsing API JSON into the expected schema.
- [ ] (Model: GPT-5.5) Review biological metadata fields and isoform handling.

## Step 11: ModelCIF, ModelPDB, And Metadata Schema Review

- [ ] (Model: GPT-5.4) Create branch
  `integration-pr-26-gpu-step-11-model-metadata`.
- [ ] (Model: GPT-5.4) Review `afdb_integration_kit/modelcif/generate.py`,
  `modelpdb/generate.py`, metadata schema changes, and template changes.
- [ ] (Model: GPT-5.4) Confirm new QA metrics are valid ModelCIF categories and
  use stable metric IDs.
- [ ] (Model: GPT-5.4) Confirm schema changes remain backward-compatible with
  existing examples and tests.
- [ ] (Model: GPT-5.4) Add tests for:
  - `cif_qa_metrics`.
  - model JSON metric enrichment.
  - struct ref clamping.
  - batch ModelCIF/ModelPDB commands if retained.
- [ ] (Model: GPT-5.5) Review ModelCIF compliance before merging.

## Step 12: Nextflow Workflow Review

- [ ] (Model: GPT-5.4) Create branch
  `integration-pr-26-gpu-step-12-nextflow`.
- [ ] (Model: GPT-5.4) Review `workflow/workflow.nf` and new multibatch
  workflows.
- [ ] (Model: GPT-5.4) Decide whether `end_to_end_with_validation_multibatch_old.nf`
  should be retained, renamed, or removed.
- [ ] (Model: GPT-5.4) Check process inputs/outputs against actual script CLI
  arguments.
- [ ] (Model: GPT-5.4) Add or update small Nextflow smoke documentation if a
  real Nextflow run is too heavy for CI.
- [ ] (Model: GPT-5.5) Review operational fit with existing pipeline usage.

## Step 13: Large Files, Examples, And Repo Shape

- [ ] (Model: GPT-5.4-mini) Create branch
  `integration-pr-26-gpu-step-13-repo-shape`.
- [ ] (Model: GPT-5.4-mini) Review deleted example zip files and newly added
  fixture/helper files.
- [ ] (Model: GPT-5.4) Decide whether binary example deletions are acceptable.
- [ ] (Model: GPT-5.4) Confirm `.gitignore` changes do not hide important
  source files.
- [ ] (Model: GPT-5.4-mini) Run `git diff --stat origin/main...HEAD` and flag
  very large additions for review.
- [ ] (Model: GPT-5.4) Commit repo-shape decisions separately from code.

## Step 14: Final Integration And Merge Readiness

- [ ] (Model: GPT-5.5) Create final review branch from the accumulated
  integration branch.
- [ ] (Model: GPT-5.5) Review all commits since `origin/main`.
- [ ] (Model: GPT-5.4-mini) Run full local verification:
  - `pytest -q`
  - `python -m compileall -q main.py afdb_integration_kit uniprot scripts tests`
  - `git diff --check origin/main...HEAD`
  - `uv lock --locked`
  - `uv export --locked --no-hashes --output-file=requirements.txt --no-dev`
  - `uv export --locked --no-hashes --extra production --output-file
    /tmp/production-req.txt --no-dev`
  - `python main.py --help`
- [ ] (Model: GPT-5.4) Build or dry-run Docker according to the final Docker
  decision.
- [ ] (Model: GPT-5.5) Write final PR review summary:
  - what was kept from PR #26.
  - what was changed for repo compatibility.
  - remaining operational risks.
  - exact verification commands and results.
- [ ] (Model: GPT-5.5) Only then merge to `main` or open the cleaned
  integration branch as a replacement PR.
