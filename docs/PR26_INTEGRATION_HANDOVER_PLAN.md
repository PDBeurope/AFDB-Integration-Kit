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

Current examples follow-up branch:

- `integration-pr-26-gpu-examples-complex-e2e`

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

## Examples Follow-Up: Complex End-to-End References

- [x] Create branch `integration-pr-26-gpu-examples-complex-e2e` from
  `integration-pr-26-gpu`.
- [x] Keep the monomer reference under `examples/colabfold_e2e/` unchanged in
  scope and add a sibling complex reference under
  `examples/colabfold_complex_e2e/`.
- [x] Select one curated homodimer and one curated heterodimer fixture that
  complete against the supplied DuckDB.
  - Homodimer: `AF-0000000065760001` -> `Q6GZX4`
  - Heterodimer: `AF-0000000300000101` -> `A0ABS2QMZ4` + `A0ABS2QMF5`
- [x] Reuse the existing stitched script flow rather than adding a Nextflow
  dependency.
- [x] Fix only the example blocker encountered in the stitched flow.
  - `afdb_integration_kit.colabfold.converter` now respects manifest-provided
    chain ranges when DuckDB metadata is used, instead of assuming every chain
    spans the full UniProt sequence length.
  - `scripts/generate_colabfold_e2e_example.py` now stages per-chain local
    residue ranges in its input manifest, which keeps complex metadata export
    consistent for the curated example set.
- [x] Add focused regression coverage for the converter and helper manifest
  behavior.
- [x] Record exact commands in
  `examples/colabfold_complex_e2e/config/commands.txt` and tool caveats in
  `examples/colabfold_complex_e2e/README.md`.
- [ ] If this branch is later merged, update the parent-branch status section
  with the merge commit and parent verification summary.

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

- [x] (Model: GPT-5.4) Create branch
  `integration-pr-26-gpu-step-6-colabfold-manifest`.
- [x] (Model: GPT-5.4) Review `afdb_integration_kit/colabfold/converter.py`
  changes.
  - Confirm gemmi parsing gives identical chain spans to the previous parser
    for representative PDBs.
    - Evidence: new synthetic multi-chain test compares
      `_chain_spans_from_pdb_gemmi()` with `_chain_spans_from_pdb_legacy()`
      and preserves the original global flattened spans.
  - Confirm PAE rounding and `orjson.OPT_SERIALIZE_NUMPY` are intended.
    - Evidence: `tests/test_colabfold_converter.py` now verifies rounded JSON
      output from `convert_file()` and confirms NumPy-backed payloads are
      written as normal JSON arrays.
  - Confirm DuckDB prefetch cache is safe across repeated conversions.
    - Evidence: partial-prefetch test now covers the previous stale-cache risk
      where a same-path DuckDB cache could be populated for only one accession
      and then incorrectly reused for a different accession later.
- [x] (Model: GPT-5.4) Review `afdb_integration_kit/manifest/resolver.py`.
  - Confirm supported model ID formats: `AF-...`, `AF_...`, and
    `AF_..._AF_...`.
    - Evidence: tests cover all three supported shapes.
  - Confirm 16-digit AF IDs are accepted if that is the production convention.
    - Evidence: classifier now explicitly accepts only 16-digit AF IDs.
  - Confirm hyphen/underscore normalization is consistent with filenames,
    manifests, and metadata.
    - Evidence: manifest-row builder test covers underscore input and
      hyphenated output `model_entity_id` rows.
  - Confirm ambiguous accession deduplication is scientifically valid.
    - Decision: when no UniProt DuckDB plus ColabFold metadata evidence is
      available, ambiguous accession sets now fail instead of choosing an
      arbitrary alphabetical accession.
- [x] (Model: GPT-5.4) Add tests for model-ID classification and manifest row
  building.
- [x] (Model: GPT-5.4) Add tests for chain mapping behavior in
  `convert_file`.
- [x] (Model: GPT-5.4) Add curated real ColabFold fixtures for later
  integration and end-to-end smoke tests.
  - Evidence: `tests/fixtures/colabfold_real_examples` now contains 3
    monomers, 3 homodimers, and 3 heterodimers with minimal raw score JSON,
    PDB, and small metadata/input files.
  - Fixture names use single AF-style `model_entity_id` values. Heterodimer
    source files that were originally named with two component AF IDs now use
    reserved single fixture IDs, while source component IDs are retained in
    `manifest.json`.
  - Chain UniProt accessions are recorded in `manifest.json` from the corrected
    merged ColabFold manifest and related corrected mapping outputs.
  - Real-fixture tests now run `convert_file()` across all 9 curated examples
    and assert output confidence lengths, PAE dimensions, global chain spans,
    source IDs, and chain accessions from `manifest.json`.
- [x] (Model: GPT-5.4-mini) Run:
  - `pytest tests/test_colabfold_converter.py -q`
    - Result: `14 passed` via
      `UV_CACHE_DIR=/tmp/uv-cache uv run pytest tests/test_colabfold_converter.py -q`.
  - new manifest resolver tests
    - Result: `4 passed` via
      `UV_CACHE_DIR=/tmp/uv-cache uv run pytest tests/test_manifest_resolver.py -q`.
  - `.venv/bin/python -m compileall -q afdb_integration_kit/colabfold
    afdb_integration_kit/manifest tests`
    - Result: passed.
  - `git diff --check`
    - Result: passed.
  - `.venv/bin/python -m pytest -q`
    - Result: `73 passed, 1 skipped, 1 warning`.
- [x] (Model: GPT-5.5) Review data-model assumptions before merging.
  - Remaining assumption: `build_colabfold_manifest()` still treats a single
    AF ID as a homodimer and emits two chains. This matches the PR branch
    implementation but is only indirectly justified by the AFCDB/ColabFold
    multimer workflow context; no broader monomer claim was added in Step 6.

Note: Step 6 has been merged back into `integration-pr-26-gpu` with merge
commit `4c4a158`. Parent-branch verification after the merge:

- `UV_CACHE_DIR=/tmp/uv-cache uv run pytest tests/test_colabfold_converter.py
  -q`: passed with `14 passed`.
- `UV_CACHE_DIR=/tmp/uv-cache uv run pytest tests/test_manifest_resolver.py
  -q`: passed with `4 passed`.
- `.venv/bin/python -m compileall -q afdb_integration_kit/colabfold
  afdb_integration_kit/manifest tests`: passed.
- `git diff --check`: passed.
- `.venv/bin/python -m pytest -q`: passed with `73 passed, 1 skipped, 1
  warning`.

The Step 7 branch was created from this verified parent after the Step 6
merge/status documentation commit. No Step 7 implementation has started yet.

## Step 7: GPU Clash/Interface Package Review

- [x] (Model: GPT-5.4) Create branch
  `integration-pr-26-gpu-step-7-gpu-analysis`.
- [x] (Model: GPT-5.4) Preserve existing toolkit behavior outside the new
  clash/interface analysis package.
  - Core toolkit imports and commands must not require Torch, fastpdb, Biotite,
    CUDA, or torch_cluster unless the user calls production analysis code.
  - Shared behavior between the current toolkit and PR #26 should remain
    functionally similar unless there is an explicit integration decision and
    test coverage for the change.
- [x] (Model: GPT-5.4) Review `afdb_integration_kit/gpu/*` as a package.
  - Confirm imports are optional and do not affect core package import.
  - Consider making `afdb_integration_kit.gpu.__init__` lazy or minimal so
    missing `fastpdb`, `torch`, or `torch_cluster` produce clear errors.
  - Treat the package as Torch-based clash/interface analysis that can run on
    CPU for small/correctness workloads and on CUDA for production-scale
    throughput.
  - Support a clear device contract for public entry points:
    `device="cpu"`, `device="cuda"`, and preferably `device="auto"`.
  - If `device="cuda"` is requested without CUDA availability, fail early with
    a clear message instead of failing later during tensor transfer.
- [x] (Model: GPT-5.4) Check licensing headers and provenance for copied GPU
  code.
- [x] (Model: GPT-5.4) Add CPU-compatible tests where possible.
  - Test schema conversion.
  - Test small clash/interface calculations if torch is available.
  - Test CPU execution explicitly for tiny synthetic `Protein` objects.
  - Test `device="auto"` resolution if added.
  - Skip CUDA-specific tests unless CUDA is available.
- [x] (Model: GPT-5.4) Check `torch_cluster` fallback behavior in
  `clashes.py`.
  - Confirm the pure-PyTorch fallback produces correct small-case results when
    `torch_cluster` is unavailable or intentionally bypassed.
  - Keep `torch_cluster` documented as an optional accelerator with
    environment-specific wheel installation.
- [x] (Model: GPT-5.5) Review performance-sensitive code for memory pressure
  and batching assumptions.
  - Document that CPU tests validate correctness and API behavior, while CUDA
    throughput and GPU memory behavior remain unverified in this CPU-only
    environment.
  - Step 7 implementation notes:
    - `afdb_integration_kit.gpu.__init__` is now lazy so plain package import
      no longer pulls `torch`, `fastpdb`, or `biotite`.
    - Public GPU analysis entry points now accept `device="auto"` and fail
      early with a clear error if `device="cuda"` is requested without CUDA.
    - `parse.py` loads `fastpdb` and Biotite lazily, so `analyze_proteins()`
      is not coupled to PDB parsing support at import time.
    - Added `tests/test_gpu_analysis.py` for lightweight import checks,
      dependency-error messaging, schema conversion, `device="auto"` resolution,
      and a CPU/fallback execution path that runs when PyTorch is installed.
    - In this sandbox `.venv` does not contain `torch`, `fastpdb`, or
      `torch_cluster`, so the CPU execution/fallback test is skip-marked and
      CUDA throughput remains unverified here.

Note: Step 7 has been merged back into `integration-pr-26-gpu` with merge
commit `0529ffa`. Parent-branch verification after the merge:

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

The Step 8 branch was created from this verified parent after the Step 7
merge/status documentation commit. Step 8 has now been merged back into the
parent and re-verified there.

## Step 8: iPSAE C++ Tool Review

- [x] (Model: GPT-5.4) Create branch
  `integration-pr-26-gpu-step-8-ipsae`.
- [x] (Model: GPT-5.4) Review `afdb_integration_kit/ipsae/ipsae_cpp.cpp` and
  `Makefile`.
  - Added direct standard-library includes used by `ipsae_cpp.cpp`
    (`<atomic>`, `<iomanip>`, `<map>`) instead of relying on transitive
    includes from other headers.
  - Updated the iPSAE `Makefile` so plain `make` prefers an existing local or
    system Eigen install (`/usr/include/eigen3`, `/usr/local/include/eigen3`,
    or a repo-local cached copy) before attempting a network fetch.
- [x] (Model: GPT-5.4) Decide whether vendored `json.hpp` is acceptable.
  - If vendored, document source/version/license.
  - If not, switch to a system/package dependency.
  - Decision: retain vendored `afdb_integration_kit/ipsae/deps/json.hpp`.
  - Rationale: it is a header-only dependency, already carries upstream SPDX
    metadata, and keeps the iPSAE build self-contained without adding a new
    package-manager dependency.
  - Documentation: added provenance/license notes to
    `afdb_integration_kit/ipsae/README.md`.
- [x] (Model: GPT-5.4) Confirm the C++ tool builds locally with available
  compilers.
  - `make check`: passed after the Makefile change, resolving Eigen from
    `/usr/include/eigen3`.
  - `make clean && make`: passed and produced `afdb_integration_kit/ipsae/ipsae_cpp`.
  - Build note: the static link produced a `libgomp.a`/`dlopen` warning from
    the host toolchain, but the build still succeeded.
- [x] (Model: GPT-5.4) Add a minimal fixture and test for expected CSV output,
  if build tooling is available.
  - Added `tests/test_ipsae_cpp.py`.
  - The test builds the tool through the iPSAE `Makefile`, runs a tiny
    two-chain batch fixture, and checks the summary CSV columns plus basic
    numeric expectations for `ipsae_AB`, `ipsae_BA`, `LIS_AB`, `LIS_BA`,
    `iptm_af`, `n0chn`, and `d0chn`.
- [x] (Model: GPT-5.4-mini) Run C++ build and smoke test where supported.
  - `.venv/bin/python -m pytest -q tests/test_ipsae_cpp.py
    tests/test_shard_analysis_metadata.py`: passed with `3 passed`.
  - `.venv/bin/python -m compileall -q afdb_integration_kit tests`: passed.
  - `git diff --check`: passed.
  - `.venv/bin/python -m pytest -q`: passed with `80 passed, 2 skipped, 1 warning`.
- [x] (Model: GPT-5.5) Review numerical correctness and threshold handling.
  - Review note: `pae_cutoff` directly gates ipSAE accumulation, while
    `dist_cutoff` only affects the reported `dist_nres1`/`dist_nres2` counts;
    it does not change the ipSAE score itself.
  - Review note: pDockQ contact detection is hard-coded at 8.0 A and LIS uses
    a hard-coded 12.0 A PAE cutoff; neither is driven by CLI thresholds.
  - Review note: missing `CB` atoms fall back to `CA` coordinates, including
    but not limited to glycine residues.
  - Review note: the batch summary schema assumes single-character chain IDs
    when constructing columns such as `ipsae_AB` and `ipsae_BA`.
  - Deferred note: `uniprot/scripts/batch_ipsae.py` still appears fork-stale
    relative to the reviewed binary path/CLI (`afdb_integration_kit/ipsae/ipsae_cpp`).
    Step 8 did not widen scope into a UniProt script rewrite because the active
    production path uses `scripts/production_pipeline.py`; revisit that script
    in Step 10 if it is still intended to be supported.

Step 8 status:

- Merged back into `integration-pr-26-gpu` with merge commit `e4cef38`.

Parent-branch verification after the merge:

- `make check` in `afdb_integration_kit/ipsae`: passed.
- `make clean && make` in `afdb_integration_kit/ipsae`: passed and rebuilt
  `ipsae_cpp`.
  - Caveat: the static link may emit a host `libgomp.a`/`dlopen` warning, but
    the build still succeeds.
- `.venv/bin/python -m pytest -q tests/test_ipsae_cpp.py
  tests/test_shard_analysis_metadata.py`: passed with `3 passed`.
- `.venv/bin/python -m compileall -q afdb_integration_kit tests`: passed.
- `git diff --check`: passed.
- `.venv/bin/python -m pytest -q`: passed with `80 passed, 2 skipped, 1 warning`.

## Step 9: Runnable ColabFold End-to-End Reference

- [x] (Model: GPT-5.4) Create branch
  `integration-pr-26-gpu-step-9-production-pipeline`.
- [x] (Model: GPT-5.4) Pivot Step 9 away from a broad production-pipeline
  review and toward a small runnable reference flow from raw ColabFold-like
  fixture outputs to final AFDB artifacts.
- [x] (Model: GPT-5.4) Use the old Nextflow workflow only as the required
  sequence/output reference; do not rewrite the Nextflow pipeline itself.
- [x] (Model: GPT-5.4) Review the active small-run path across
  `scripts/prepare_inputs.py`, `scripts/production_pipeline.py`,
  `afdb_integration_kit.colabfold.converter`, `uniprot/scripts/*`, and the
  example fixtures.
  - Outcome: keep `scripts/production_pipeline.py` unchanged for now; it is
    broader than needed for the Step 9 proof and not the cleanest reproduction
    path for a tiny local example.
  - Outcome: add a narrow repo-owned helper,
    `scripts/generate_colabfold_e2e_example.py`, that stages a few curated
    fixture inputs and runs the existing toolkit commands in the Nextflow
    order.
- [x] (Model: GPT-5.4) Confirm the example path works entirely offline apart
  from access to the supplied local DuckDB.
  - No network calls are required.
  - The example uses three curated monomer fixtures:
    `AF-0000000300000001` (`O00400`),
    `AF-0000000300000002` (`O64637`),
    `AF-0000000300000003` (`Q9TVL3`).
- [x] (Model: GPT-5.4) Fix only direct blockers in the selected flow.
  - Fixed `uniprot/scripts/export_modelcif_input.py` by restoring the missing
    `import os`.
  - Preserved Step 1-8 behavior; no broad UniProt/API rewrite was started.
- [x] (Model: GPT-5.4) Populate `examples/colabfold_e2e/` with generated
  reference outputs and exact commands.
  - Added `examples/colabfold_e2e/README.md`.
  - Added `examples/colabfold_e2e/config/commands.txt` and
    `examples/colabfold_e2e/run_summary.json`.
  - Added normalized staged inputs, converted AFDB JSONs, merged manifests,
    model/chain metadata JSONs and batch files, ModelCIF input JSONs, mmCIFs,
    DSSP outputs, enriched PDBs, and BCIF outputs.
- [x] (Model: GPT-5.4) Record optional-tool behavior clearly.
  - `mkdssp` was available locally and used for DSSP.
  - Mol* `cif2bcif` was not on `PATH`.
  - The explicit Biotite backend worked for plain ModelCIF files but failed on
    the DSSP-enriched CIFs; the committed example BCIF files were therefore
    generated from the pre-DSSP ModelCIF files instead, and that fallback is
    recorded in `run_summary.json`.
- [x] (Model: GPT-5.4) Add focused tests for the helper/manifests rather than a
  heavy full end-to-end pytest run against the large DuckDB.
  - Added `tests/test_generate_colabfold_e2e_example.py`.
- [x] (Model: GPT-5.4-mini) Verify Step 9.
  - `.venv/bin/python main.py --help`: passed.
  - `.venv/bin/python scripts/production_pipeline.py --help`: passed.
  - `.venv/bin/python -m afdb_integration_kit.colabfold.converter --help`:
    passed.
  - `.venv/bin/python uniprot/scripts/export_model_metadata.py --help`:
    passed.
  - `.venv/bin/python uniprot/scripts/export_chain_metadata.py --help`:
    passed.
  - `.venv/bin/python uniprot/scripts/export_modelcif_input.py --help`:
    passed.
  - `.venv/bin/python uniprot/scripts/combine_metadata.py --help`: passed.
  - `.venv/bin/python scripts/generate_colabfold_e2e_example.py --duckdb
    /mnt/disks/toolkit-data/uniprot_extract_2025_04_merged_5way/db/uniprot_2025_04_merged_5way.duckdb
    --output-dir examples/colabfold_e2e`: passed.
  - Merge commit on the parent branch: `ff96f0a` (`merge
    integration-pr-26-gpu-step-9-production-pipeline`).

## Examples Follow-Up Track: Complex End-To-End References

The broad PR #26 step checklist below is deferred for now. The current priority
is to use this repository as a toolkit and script collection, not to support
every possible orchestration style or fully certify the old Nextflow pipeline.

Step 9 established a small monomer reference at `examples/colabfold_e2e/`.
Future work should keep that reference reproducible and then add comparable
runnable references for curated ColabFold-style complexes.

- [ ] (Model: GPT-5.4) Create branch
  `integration-pr-26-gpu-examples-complex-e2e` from `integration-pr-26-gpu`.
- [ ] (Model: GPT-5.4) Re-run or inspect the committed monomer reference under
  `examples/colabfold_e2e/` and confirm the README, commands, run summary, and
  generated artifacts are coherent.
- [ ] (Model: GPT-5.4) Use the old Nextflow workflow only as the required
  sequence/output reference.
  - Do not rewrite the Nextflow pipeline.
  - Do not require users to use Nextflow.
  - Treat the e2e examples as script-stitching proofs that individual toolkit
    commands can be composed.
- [ ] (Model: GPT-5.4) Add a small homodimer example using the curated fixtures
  under `tests/fixtures/colabfold_real_examples/homodimers`.
- [ ] (Model: GPT-5.4) Add a small heterodimer example using the curated
  fixtures under `tests/fixtures/colabfold_real_examples/heterodimers`.
- [ ] (Model: GPT-5.4) Use the supplied local DuckDB for metadata exports:
  `/mnt/disks/toolkit-data/uniprot_extract_2025_04_merged_5way/db/uniprot_2025_04_merged_5way.duckdb`.
- [ ] (Model: GPT-5.4) Fix only blockers that prevent the stitched script
  sequence from producing reference artifacts.
  - Keep broad UniProt/API/Nextflow/ModelCIF reviews deferred unless they
    directly block the examples.
- [ ] (Model: GPT-5.4-mini) Verify with focused helper tests and at least one
  manual/example generation command.
  - Record optional external-tool behavior clearly, including DSSP, Mol*
    `cif2bcif`, Biotite fallback, iPSAE, Torch, and DuckDB requirements.

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
