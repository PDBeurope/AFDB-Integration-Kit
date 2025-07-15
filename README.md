# AFDB Integration Toolkit

A comprehensive toolkit for integrating structural models into the AlphaFold Database (AFDB). This toolkit provides essential tools and workflows to prepare, validate, and format molecular structure data for seamless integration with AFDB infrastructure.

## Table of Contents

- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Docker Usage](#docker-usage)
- [Nextflow Workflow](#nextflow-workflow)
- [File Structure Requirements](#file-structure-requirements)
- [Contributing](#contributing)
- [License](#license)
- [Support](#support)

## Features

- **ModelCIF Generation**: Convert PDB files to mmCIF format with metadata integration
- **Binary CIF Conversion**: Efficient conversion from mmCIF to Binary CIF (BCIF) format
- **Secondary Structure Assignment**: DSSP-based secondary structure annotation
- **Automated Workflows**: Nextflow-based end-to-end processing pipelines
- **Docker Support**: Containerized execution for reproducible results
- **Validation Tools**: Built-in testing and validation utilities

## Prerequisites

- Python 3.12+
- Node.js 18+ (for Mol* CLI)
- Docker (optional, for containerized execution)
- Nextflow (optional, for workflow automation)

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/PDBeurope/AFDB-Integration-Kit
cd AFDB-Integration-Kit
```

### 2. Install UV (Python Package Manager)

UV is used to manage Python dependencies and virtual environments.

**macOS/Linux:**
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**Windows:**
```bash
powershell -c "irm https://astral.sh/uv/install.ps1 | iex"
```

**Alternative installation methods:**
```bash
# Using pip
pip install uv

# Using conda
conda install -c conda-forge uv
```

### 3. Install Mol* CLI

If you use nvm (Node Version Manager):
```bash
nvm use  # Uses the version specified in .nvmrc
npm install -g molstar
```

Without nvm:
```bash
npm install -g molstar
```

### 4. Install DSSP

We use the modern DSSP implementation by the PDB-REDO team:

```bash
# Clone and build DSSP
git clone https://github.com/PDB-REDO/dssp.git
cd dssp
mkdir build
cd build
cmake ..
make
sudo make install
```

For detailed installation instructions, visit: https://github.com/PDB-REDO/dssp

### 5. Install Gemmi

Gemmi is a structural biology library with command-line tools:

**Using conda:**
```bash
conda install -c conda-forge gemmi
```

**Using pip:**
```bash
pip install gemmi
```

**From source:**
```bash
git clone https://github.com/project-gemmi/gemmi.git
cd gemmi
make
```

For more installation options, see: https://gemmi.readthedocs.io/en/latest/install.html

### 6. Download mmCIF Dictionary (Required for ModelCIF Generator)

The ModelCIF tool requires the updated mmCIF dictionary for validation. Download it to your project directory:

```bash
# Download the mmCIF dictionary
curl -o mmcif_ma.dic https://raw.githubusercontent.com/ihmwg/ModelCIF/refs/heads/master/dist/mmcif_ma.dic
```

**Note:** This step is automatically handled in the Docker environment, but is required for local installations.

### 7. Install Nextflow (Optional)

For workflow automation:

```bash
# Using curl
curl -s https://get.nextflow.io | bash

# Make executable and add to PATH
chmod +x nextflow
sudo mv nextflow /usr/local/bin/
```

### 8. Install Docker (Optional)

For containerized execution:
- **macOS/Windows**: Download Docker Desktop from https://www.docker.com/products/docker-desktop
- **Linux**: Follow instructions at https://docs.docker.com/engine/install/

## Quick Start

### Verify Installation

Test that all dependencies are correctly installed:

```bash
uv run main.py test
```

This command will validate your environment and report any missing dependencies.

### Basic Usage Example

```bash
# Generate ModelCIF
uv run main.py run-modelcif-gen \
    -p input/AF-0000000000000001-model-v1.pdb \
    -m input/AF-0000000000000001-v1.cif.json \
    -o output/AF-0000000000000001-model-v1.cif

# Convert to BCIF
uv run main.py run-cif2bcif \
    -i input/AF-0000000000000001-model-v1.cif \
    -o output/AF-0000000000000001-model-v1.bcif

# Add secondary structure annotation
uv run main.py run-dssp \
    -i input/AF-0000000000000001-model-v1.cif \
    -o output/AF-0000000000000001-model-v1.cif
```

## Usage

### ModelCIF Generator

Converts PDB files to mmCIF format with integrated metadata.

**Requirements:**
- Input PDB file
- Metadata JSON file conforming to the schema: `afdb_integration_kit/modelcif/resources/schema.json`
- mmCIF dictionary file (mmcif_ma.dic) in the project directory

**Important:** For local installations, ensure you have downloaded the mmCIF dictionary:
```bash
curl -o mmcif_ma.dic https://raw.githubusercontent.com/ihmwg/ModelCIF/refs/heads/master/dist/mmcif_ma.dic
```

**Command:**
```bash
uv run main.py run-modelcif-gen -p <pdb_file> -m <metadata_json> -o <output_cif>
```

**Parameters:**
- `-p, --pdb`: Input PDB file path
- `-m, --metadata`: Metadata JSON file path
- `-o, --output`: Output mmCIF file path

### CIF to BCIF Converter

Converts mmCIF files to Binary CIF format for efficient storage and transmission.

**Command:**
```bash
uv run main.py run-cif2bcif -i <input_cif> -o <output_bcif>
```

**Parameters:**
- `-i, --input`: Input mmCIF file path
- `-o, --output`: Output BCIF file path

### DSSP Secondary Structure Assignment

Assigns secondary structure annotations based on atomic coordinates.

**Command:**
```bash
uv run main.py run-dssp -i <input_cif> -o <output_cif>
```

**Parameters:**
- `-i, --input`: Input mmCIF file path
- `-o, --output`: Output annotated mmCIF file path

## Docker Usage

### Build Docker Image

```bash
docker build -t afdb-kit .
```

### Run Tools in Docker

**ModelCIF Generator:**
```bash
docker run \
    -v "$PWD/input:/input" \
    -v "$PWD/output:/output" \
    -w /workspace \
    -v "$PWD:/workspace" \
    afdb-kit uv run main.py run-modelcif-gen \
        -p /input/AF-0000000000000001-model-v1.pdb \
        -m /input/AF-0000000000000001-v1.cif.json \
        -o /output/AF-0000000000000001-model-v1.cif
```

**CIF to BCIF Converter:**
```bash
docker run \
    -v "$PWD/input:/input" \
    -v "$PWD/output:/output" \
    -w /workspace \
    -v "$PWD:/workspace" \
    afdb-kit uv run main.py run-cif2bcif \
        -i /input/AF-0000000000000001-model-v1.cif \
        -o /output/AF-0000000000000001-model-v1.bcif
```

**DSSP Processing:**
```bash
docker run \
    -v "$PWD/input:/input" \
    -v "$PWD/output:/output" \
    -w /workspace \
    -v "$PWD:/workspace" \
    afdb-kit uv run main.py run-dssp \
        -i /input/AF-0000000000000001-model-v1.cif \
        -o /output/AF-0000000000000001-model-v1.cif
```

## Nextflow Workflow

### End-to-End Processing

Run the complete workflow using the provided script:

```bash
./run_workflow.sh
```

### Input Requirements

The Nextflow workflow requires an input list file at `input/input.txt` containing the entries to process. Each entry should be on a new line:

```
AF-0001234567890123
AF-0001234567890124
AF-0001234567890125
AF-0001234567890126
```

**Example input.txt:**
```bash
# Create the input list file
cat > input/input.txt << EOF
AF-0001234567890123
AF-0001234567890124
AF-0001234567890125
EOF
```

### Workflow Features

- **Resumable**: Uses `-resume` flag to continue from previous checkpoints
- **Cached**: Maintains state in `.nextflow` directory
- **Dependency Management**: Automatically handles tool dependencies
- **Parallel Processing**: Processes multiple files concurrently

### Important Notes

- Mount the `.nextflow` directory to preserve workflow state
- Ensure proper input/output directory mounting
- The workflow runs in resume mode by default

## File Structure Requirements

### Input Directory Structure

The toolkit expects files to be organized in a specific hierarchical structure:

```
input/
├── 0001/
│   ├── 2345/
│   │   ├── 6789/
│   │   │   ├── 0123/
│   │   │   │   ├── AF-0001234567890123-model-v1.pdb
│   │   │   │   └── AF-0001234567890123-v1.cif.json
```

### Directory Structure Rules

1. **Extract 16-digit numeric ID**: From `AF-0001234567890123-model-v1.pdb` → `0001234567890123`
2. **Split into 4-digit segments**: `0001`, `2345`, `6789`, `0123`
3. **Create nested directories**: `0001/2345/6789/0123/`
4. **Place files in final directory**: Both PDB and JSON files

### Output Structure

The workflow automatically creates corresponding output directories following the same structure:

```
output/
├── 0001/
│   ├── 2345/
│   │   ├── 6789/
│   │   │   ├── 0123/
│   │   │   │   ├── AF-0001234567890123-model-v1.cif
│   │   │   │   └── AF-0001234567890123-model-v1.bcif
```

## Troubleshooting

### Common Issues

1. **Missing Dependencies**: Run `uv run main.py test` to identify missing components
2. **Permission Errors**: Ensure Docker has proper access to mounted directories
3. **File Not Found**: Verify input files follow the required directory structure
4. **Memory Issues**: For large datasets, consider adjusting Docker memory limits
5. **ModelCIF Validation Errors**: Ensure `mmcif_ma.dic` is present in the project directory (automatically handled in Docker)
6. **Nextflow Workflow Errors**: Ensure `input/input.txt` exists and contains valid entry IDs

### Getting Help

- Check the [Issues](https://github.com/PDBeurope/AFDB-Integration-Kit/issues) page
- Validate your metadata JSON against the provided schema


## License

This project is licensed under the [CC0 1.0 Universal](LICENSE) - see the LICENSE file for details.

## Support

For support and questions:

- **Issues**: [GitHub Issues](https://github.com/PDBeurope/AFDB-Integration-Kit/issues)
- **Email**: afdbhelp@ebi.ac.uk


---
