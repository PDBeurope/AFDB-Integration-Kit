from __future__ import annotations

OLIGOMERIC_STATES = [
    "monomer",
    "dimer",
    "trimer",
    "tetramer",
    "pentamer",
    "hexamer",
    "heptamer",
    "octamer",
    "nonamer",
    "decamer",
    "oligomer",
]

ASSEMBLY_TYPES = ["Homo", "Hetero"]

CHAIN_ENTRY_SCHEMA = {
    "type": "object",
    "additionalProperties": False,
    "required": [
        "uniqueId",
        "toolUsed",
        "modelCreatedDate",
        "modelEntityId",
        "isComplex",
        "complexName",
        "uniprotAccession",
        "uniprotId",
        "uniprotDescription",
        "geneSynonyms",
        "taxId",
        "organismScientificName",
        "organismCommonNames",
        "organismSynonyms",
        "assemblyType",
        "oligomericState",
        "sequence",
        "sequenceChecksum",
        "sequenceVersionDate",
        "sequenceStart",
        "sequenceEnd",
        "isIsoform",
        "isFragment",
        "isUniProt",
        "isUniProtReferenceProteome",
        "isUniProtReviewed",
        "globalMetricValue",
        "fractionPlddtVeryLow",
        "fractionPlddtLow",
        "fractionPlddtConfident",
        "fractionPlddtVeryHigh",
        "latestVersion",
        "allVersions",
        "providerId",
        "entityType",
        "isAMdata",
    ],
    "properties": {
        "uniqueId": {"type": "string", "pattern": r"^AF-\d{16}_v[1-9]\d*_[A-Za-z0-9]+$"},
        "toolUsed": {"type": "string", "minLength": 1},
        "modelCreatedDate": {"type": "string", "format": "date-time"},
        "modelEntityId": {"type": "string", "pattern": r"^AF-\d{16}-v[1-9]\d*$"},
        "isComplex": {"type": "boolean"},
        "complexName": {"type": ["string", "null"]},
        "complexComposition": {
            "type": "string",
            "pattern": r"^[A-Za-z0-9-]+_[1-9]\d*(,[A-Za-z0-9-]+_[1-9]\d*)*$",
        },
        "complexPredictionAccuracy_ipTM": {"type": "number", "minimum": 0.0, "maximum": 1.0},
        "uniprotAccession": {"type": "string", "minLength": 1},
        "uniprotId": {"type": "string", "minLength": 1},
        "uniprotDescription": {"type": "string", "minLength": 1},
        "gene": {"type": "string", "minLength": 1, "pattern": r".*\S.*"},
        "geneSynonyms": {"type": "array", "items": {"type": "string"}},
        "taxId": {"type": "integer", "minimum": 1},
        "organismScientificName": {"type": "string", "minLength": 1},
        "organismCommonNames": {"type": "array", "items": {"type": "string"}},
        "organismSynonyms": {"type": "array", "items": {"type": "string"}},
        "assemblyType": {"type": ["string", "null"], "enum": [None, *ASSEMBLY_TYPES]},
        "oligomericState": {"type": ["string", "null"], "enum": [None, *OLIGOMERIC_STATES]},
        "sequence": {"type": "string", "minLength": 1, "pattern": r"^[A-Z]+$"},
        "sequenceChecksum": {"type": "string", "pattern": r"^[a-fA-F0-9]{32}$"},
        "sequenceVersionDate": {"type": "string", "format": "date-time"},
        "sequenceStart": {"type": "integer", "minimum": 1},
        "sequenceEnd": {"type": "integer", "minimum": 1},
        "isIsoform": {"type": "boolean"},
        "isFragment": {"type": "boolean"},
        "isUniProt": {"type": "boolean"},
        "isUniProtReferenceProteome": {"type": "boolean"},
        "isUniProtReviewed": {"type": "boolean"},
        "globalMetricValue": {"type": "number", "minimum": 0.0, "maximum": 100.0},
        "fractionPlddtVeryLow": {"type": "number", "minimum": 0.0, "maximum": 1.0},
        "fractionPlddtLow": {"type": "number", "minimum": 0.0, "maximum": 1.0},
        "fractionPlddtConfident": {"type": "number", "minimum": 0.0, "maximum": 1.0},
        "fractionPlddtVeryHigh": {"type": "number", "minimum": 0.0, "maximum": 1.0},
        "latestVersion": {"type": "integer", "minimum": 1},
        "allVersions": {"type": "array", "minItems": 1, "items": {"type": "integer", "minimum": 1}},
        "providerId": {"type": "string", "minLength": 1},
        "entityType": {"type": "string", "minLength": 1},
        "isAMdata": {"type": "boolean"},
    },
    "allOf": [
        {
            "if": {"properties": {"isComplex": {"const": True}}},
            "then": {
                "required": [
                    "complexComposition",
                    "complexPredictionAccuracy_ipTM",
                ],
                "properties": {
                    "assemblyType": {"type": "string", "enum": ASSEMBLY_TYPES},
                    "oligomericState": {"type": "string", "enum": OLIGOMERIC_STATES},
                },
            },
            "else": {
                "properties": {
                    "assemblyType": {"type": "null"},
                    "oligomericState": {"type": "null"},
                },
                "not": {
                    "anyOf": [
                        {"required": ["complexComposition"]},
                        {"required": ["complexPredictionAccuracy_ipTM"]},
                    ]
                },
            },
        }
    ],
}

MODEL_ENTRY_SCHEMA = {
    "type": "object",
    "additionalProperties": False,
    "required": [
        "modelEntityId",
        "latestVersion",
        "providerId",
        "isComplex",
        "uniprotAccession",
        "uniprotDescription",
        "isUniProtReferenceProteome",
        "isUniProtReviewed",
        "isIsoform",
        "organismScientificName",
        "taxId",
        "globalMetricValue",
        "isAMdata",
    ],
    "properties": {
        "modelEntityId": {"type": "string", "pattern": r"^AF-\d{16}$"},
        "latestVersion": {"type": "integer", "minimum": 1},
        "providerId": {"type": "string", "minLength": 1},
        "isComplex": {"type": "boolean"},
        "complexName": {"type": "string", "minLength": 1},
        "assemblyType": {"type": "string", "enum": ASSEMBLY_TYPES},
        "oligomericState": {"type": "string", "enum": OLIGOMERIC_STATES},
        "complexPredictionAccuracy_ipTM": {"type": "number", "minimum": 0.0, "maximum": 1.0},
        "uniprotAccession": {"type": "array", "minItems": 1, "items": {"type": "string", "minLength": 1}},
        "uniprotDescription": {"type": "array", "minItems": 1, "items": {"type": "string", "minLength": 1}},
        "isUniProtReferenceProteome": {"type": "boolean"},
        "isUniProtReviewed": {"type": "boolean"},
        "isIsoform": {"type": "boolean"},
        "organismScientificName": {"type": "array", "minItems": 1, "items": {"type": "string", "minLength": 1}},
        "gene": {
            "type": "array",
            "minItems": 1,
            "items": {"type": "string", "minLength": 1, "pattern": r".*\S.*"},
        },
        "taxId": {"type": "array", "minItems": 1, "items": {"type": "integer", "minimum": 1}},
        "sequenceStart": {"type": "integer", "minimum": 1},
        "sequenceEnd": {"type": "integer", "minimum": 1},
        "globalMetricValue": {"type": "number", "minimum": 0.0, "maximum": 100.0},
        "isAMdata": {"type": "boolean"},
    },
    "allOf": [
        {
            "if": {"properties": {"isComplex": {"const": True}}},
            "then": {
                "required": [
                    "assemblyType",
                    "oligomericState",
                    "complexPredictionAccuracy_ipTM",
                ],
            },
            "else": {
                "not": {
                    "anyOf": [
                        {"required": ["assemblyType"]},
                        {"required": ["oligomericState"]},
                        {"required": ["complexPredictionAccuracy_ipTM"]},
                        {"required": ["complexName"]},
                    ]
                },
            },
        },
        {
            "if": {"required": ["sequenceStart"]},
            "then": {"required": ["sequenceEnd"]},
        },
        {
            "if": {"required": ["sequenceEnd"]},
            "then": {"required": ["sequenceStart"]},
        },
    ],
}

METADATA_BATCH_SCHEMAS = {
    "chain": CHAIN_ENTRY_SCHEMA,
    "model": MODEL_ENTRY_SCHEMA,
}
