# AlphaFold Interface Annotation Schema

This folder packages the JSON schema, examples, and tooling we need to describe AlphaFold Complex Database residue–residue interfaces. It adapts the FunPDBe schema so it fits AF IDs (`AF-################`) while keeping the residue/site layout that PDBe-KB expects.

## Files
- `interface_schema.json` – JSON schema defining the entry metadata (`af_id`, release info), per-chain residue records (`residue_number`, `aa_type`, `site_data`), and site definitions. `af_id` must match the AlphaFold format `AF-################`.
- `interface_schema_example.json` – minimal single-residue example showing the required fields.
- `interface_mock_data.json` – concrete three-pair mock dataset: two A/B contacts and one C/D contact, with interaction details grouped inside each site's `additional_site_annotations.interactions` block.
- `validate_interface_json.py` – helper script that runs `jsonschema.validate` to make sure a JSON file matches the schema.

## Data Model Highlights
- Each AlphaFold complex entry contains a single `af_id`, resource metadata, and `chains`.
- Every chain lists only the residues that participate in interactions. Each residue has a `residue_number` (matching the UniProt/AlphaFold index), three-letter `aa_type`, and one or more `site_data` records. `site_data.site_id_ref` links the residue back to a site. `site_data` objects must contain `confidence_classification` (enum: `high`, `medium`, `low`, `null`, `curated`). Even for binary contacts you can pick a single value (e.g. `high`). In the mock data `raw_score` carries the inter-residue distance and `raw_score_unit` is `A`, but you can replace these with any residue-level metric/units your pipeline produces.
- `sites` describe interaction patches (interfaces). The mock data stores residue-pair details inside `additional_site_annotations.interactions`, so pair metadata exists once per site while the residue lists remain chain-centric. Each entry in `interactions` uses the keys `chain_1`, `res_1_label`, `aa_1_type`, `chain_2`, `res_2_label`, `aa_2_type`, and `distance_angstrom`; extend that object with more fields if you capture additional pair-level metrics.

## Validating JSON
1. Install dependencies (once): `python3 -m pip install jsonschema`.
2. Run the validator from this directory:
   ```bash
   python3 validate_interface_json.py interface_schema.json interface_mock_data.json
   ```
   Replace `interface_mock_data.json` with any file you want to check. A success message confirms the JSON satisfies the schema.

## Workflow
1. Duplicate `interface_mock_data.json` as a starting template for a new complex.
2. Adjust entry-level metadata (`af_id`, release date, links, etc.).
3. For each interface, add residues under their chains and list the residue pairs once inside the matching site's `additional_site_annotations.interactions` array.
4. Validate with the script before committing or sharing.
