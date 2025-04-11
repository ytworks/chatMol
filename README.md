# ChatMol - Molecular Property Calculator MCP Server

A Model Context Protocol (MCP) server for calculating various molecular properties from chemical structures (SMILES notation). It processes data in CSV format, calculates properties such as molecular weight, lipophilicity, hydrogen bonding capability, molecular formula, and adds the results. When used with Claude Desktop, it enables simple analysis of chemical data.

## Overview

This MCP server calculates various molecular properties from SMILES notations in CSV data and adds them as new columns:

- **Basic Properties**: Molecular weight, exact molecular weight, heavy atom molecular weight, molecular formula
- **Hydrophilicity/Hydrophobicity**: LogP (octanol/water partition coefficient), molar refractivity, topological polar surface area (TPSA)
- **Hydrogen Bonding**: Number of hydrogen bond donors, number of hydrogen bond acceptors
- **Atom & Bond Information**: Heavy atom count, heteroatom count, rotatable bond count
- **Ring Structure Information**: Number of aromatic rings, aliphatic rings, saturated rings, etc.
- **Complexity Indicators**: Number of stereogenic centers, spiro atoms, bridgehead atoms, etc.
- **Drug-likeness Filters**: Lipinski's Rule of Five, Veber's Rules, PAINS filter, etc.

## Requirements

- Python 3.10 or higher
- RDKit 2024.9.6
- pandas 2.2.3
- mcp (Model Context Protocol) package v1.2.0 or higher
- Claude Desktop
- uv (Python package manager, optional)

## Installation

### 1. Install Dependencies

```bash
# Using uv
uv pip install "rdkit==2024.9.6" "pandas==2.2.3" "mcp[cli,server]>=1.2.0"

# Using pip
pip install "rdkit==2024.9.6" "pandas==2.2.3" "mcp[cli,server]>=1.2.0"
```

## Running Tests

To verify that ChatMol is working properly, you can run the test suite:

```bash
# Run all tests
python -m pytest

# Run specific test modules
python -m pytest tests/test_properties.py
python -m pytest tests/test_io.py

# Run tests with verbose output
python -m pytest -v

# Run a specific test function
python -m pytest tests/test_properties.py::TestMolecularProperties::test_molecular_weight_calculation
```

### Test Coverage

The test suite includes:
- Validation of molecular property calculations against known values
- Verification that all descriptors can be calculated
- Proper handling of invalid input
- Performance tracking of property calculations across diverse molecular structures

## Complete List of Calculable Molecular Descriptors

ChatMol can calculate the following molecular descriptors:

### Basic Properties
| Property            | Description                                                                 |
| ------------------- | --------------------------------------------------------------------------- |
| `molecular_weight`  | Average molecular weight based on average atomic masses of elements         |
| `exact_mol_wt`      | Exact molecular weight considering isotopic composition (monoisotopic mass) |
| `heavy_atom_mol_wt` | Molecular weight ignoring hydrogens                                         |
| `formula`           | Chemical formula of the molecule                                            |

### Lipophilicity/Hydrophilicity
| Property | Description                                                                       |
| -------- | --------------------------------------------------------------------------------- |
| `logp`   | Partition coefficient LogP (1-octanol/water), estimated by Wildman-Crippen method |
| `mol_mr` | Molar refractivity, related to molecular polarizability                           |
| `tpsa`   | Topological Polar Surface Area, sum of polar atom surface areas                   |

### Surface Properties
| Property     | Description                       |
| ------------ | --------------------------------- |
| `labute_asa` | Labute's Approximate Surface Area |

### Hydrogen Bonding and Atom Counts
| Property                | Description                          |
| ----------------------- | ------------------------------------ |
| `num_h_donors`          | Number of hydrogen bond donors       |
| `num_h_acceptors`       | Number of hydrogen bond acceptors    |
| `num_rotatable_bonds`   | Number of rotatable bonds            |
| `heavy_atom_count`      | Number of heavy (non-hydrogen) atoms |
| `num_hetero_atoms`      | Number of heteroatoms (non-carbon)   |
| `no_count`              | Number of nitrogen and oxygen atoms  |
| `nhoh_count`            | Number of NH and OH groups           |
| `num_valence_electrons` | Total number of valence electrons    |

### Ring Information
| Property                     | Description                                               |
| ---------------------------- | --------------------------------------------------------- |
| `num_aromatic_rings`         | Number of aromatic rings                                  |
| `num_aliphatic_rings`        | Number of aliphatic rings                                 |
| `num_saturated_rings`        | Number of saturated rings                                 |
| `num_aromatic_carbocycles`   | Number of aromatic rings where all atoms are carbon       |
| `num_aromatic_heterocycles`  | Number of aromatic rings containing heteroatoms           |
| `num_aliphatic_carbocycles`  | Number of aliphatic rings consisting only of carbon atoms |
| `num_aliphatic_heterocycles` | Number of aliphatic rings containing heteroatoms          |
| `num_saturated_carbocycles`  | Number of saturated carbon rings                          |
| `num_saturated_heterocycles` | Number of saturated heterocyclic rings                    |
| `ring_count`                 | Total number of ring structures                           |

### Bond Information
| Property        | Description                                   |
| --------------- | --------------------------------------------- |
| `fraction_csp3` | Fraction of carbon atoms in sp³ hybridization |

### Graph Indices
| Property          | Description                                                     |
| ----------------- | --------------------------------------------------------------- |
| `balaban_j`       | Balaban's molecular connectivity index J                        |
| `bertz_ct`        | Bertz complexity index for molecular structure                  |
| `ipc`             | Information content descriptor of the molecular graph           |
| `hall_kier_alpha` | Hall-Kier alpha parameter for molecular correction              |
| `kappa1`          | Kappa shape index 1 (degree of molecular branching)             |
| `kappa2`          | Kappa shape index 2 (spatial extent, planar)                    |
| `kappa3`          | Kappa shape index 3 (spatial extent, three-dimensional)         |
| `chi0`            | Molecular connectivity index (zeroth-order)                     |
| `chi1`            | Molecular connectivity index (first-order)                      |
| `chi0v`           | Molecular connectivity index considering valence (zeroth-order) |
| `chi1v`           | Molecular connectivity index considering valence (first-order)  |

### Drug-likeness
| Property | Description                                               |
| -------- | --------------------------------------------------------- |
| `qed`    | Quantitative Estimation of Drug-likeness (score from 0-1) |

### Drug-likeness Filters
| Filter               | Description                                                                          |
| -------------------- | ------------------------------------------------------------------------------------ |
| `lipinski_pass`      | Lipinski's Rule of Five (MW≤500, LogP≤5, HBD≤5, HBA≤10)                              |
| `veber_pass`         | Veber's Rules (TPSA≤140 Å², RotBonds≤10)                                             |
| `ghose_pass`         | Ghose Filter (160≤MW≤480, -0.4≤LogP≤5.6, 20≤atoms≤70, 40≤MR≤130)                     |
| `egan_pass`          | Egan Filter (LogP≤5.88, TPSA≤131.6)                                                  |
| `muegge_pass`        | Muegge Filter (200≤MW≤600, -2≤LogP≤5, TPSA≤150, rings≤7, HBA≤10, HBD≤5, RotBonds<15) |
| `pains_free`         | PAINS filter (screens for pan-assay interference compounds)                          |
| `all_filters_passed` | Compound passes all drug-likeness filters                                            |

### Fragment Analysis
ChatMol also calculates numerous fragment-based descriptors that count the occurrences of specific functional groups within a molecule. These include:

- `fr_Al_COO` - Aliphatic carboxylic acid
- `fr_Al_OH` - Aliphatic hydroxyl
- `fr_Al_OH_noTert` - Aliphatic hydroxyl excluding tertiary
- `fr_ArN` - Aromatic nitrogen
- `fr_Ar_COO` - Aromatic carboxylic acid
- `fr_Ar_N` - Aromatic nitrogen
- `fr_Ar_NH` - Aromatic amine
- `fr_Ar_OH` - Aromatic hydroxyl (phenol)
- `fr_COO` - Carboxylic acid
- `fr_COO2` - Carboxylic acid derivative
- `fr_C_O` - Carbonyl group
- `fr_C_O_noCOO` - Carbonyl excluding carboxylic acids
- `fr_C_S` - Carbon-sulfur bonds
- `fr_benzene` - Benzene rings
- `fr_ester` - Ester groups
- `fr_ether` - Ether groups
- `fr_phenol` - Phenol groups
- `fr_ketone` - Ketone groups
- `fr_nitro` - Nitro groups

*Note: Over 70 fragment descriptors are available. This is a subset of the most commonly used ones.*

## Usage

### 1. Using with Claude Desktop

#### Edit the configuration file

Edit the Claude Desktop configuration file to add this MCP server. The configuration file path is:

- macOS: `~/Library/Application Support/Claude/claude_desktop_config.json`
- Windows: `%APPDATA%\Claude\claude_desktop_config.json`

Add the following JSON configuration (within the existing `mcpServers` object):

```json
{
  "mcpServers": {
    "ChatMol Molecular Property Calculator": {
      "command": "uv",
      "args": [
        "run",
        "--with",
        "mcp[cli]",
        "--with",
        "rdkit",
        "--with",
        "pandas",
        "mcp",
        "run",
        "/path/to/chatMol/server.py"
      ]
    }
  }
}
```

Notes:
- If the `uv` command is not in your environment path, use an absolute path (e.g.: `/path/to/uv`).
- Replace `/path/to/chatMol/server.py` with the absolute path to this script.
- Always use absolute paths, not relative paths.

### 2. Available Tools

This MCP server provides the following tools:

#### calculate_molecular_properties

- Description: Calculates molecular properties from SMILES notation
- Parameters:
  - `input_data`: SMILES string or CSV data to process (required)
  - `input_type`: Type of input data - "smiles" (single SMILES string) or "csv" (CSV data)
  - `smiles_column`: Column name containing SMILES structures (if omitted, uses the rightmost column)

#### get_available_features

- Description: Returns a list of all molecular properties that can be calculated
- Parameters: none

### Examples

In Claude Desktop, you can use it as follows:

```
Please calculate the molecular properties for this CSV data:

ID,Name,SMILES
1,Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O
2,Paracetamol,CC(=O)NC1=CC=C(C=C1)O
3,Ibuprofen,CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
```

To calculate multiple specific properties simultaneously:

```
Please calculate molecular weight, LogP, hydrogen bond donors, hydrogen bond acceptors, and molecular formula from this CSV data:

ID,Name,SMILES
1,Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O
2,Paracetamol,CC(=O)NC1=CC=C(C=C1)O
3,Ibuprofen,CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
```

To specify a column name:

```
Please calculate molecular properties from the smiles_col column in this CSV data:

ID,Name,smiles_col
1,Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O
2,Paracetamol,CC(=O)NC1=CC=C(C=C1)O
3,Ibuprofen,CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
```

To get a list of available properties:

```
What molecular properties can you calculate?
```

### Examples of Available Molecular Properties

- `molecular_weight`: Molecular weight
- `logp`: Octanol/water partition coefficient
- `num_h_donors`: Number of hydrogen bond donors
- `num_h_acceptors`: Number of hydrogen bond acceptors
- `formula`: Molecular formula
- `tpsa`: Topological polar surface area
- `num_rotatable_bonds`: Number of rotatable bonds
- `num_aromatic_rings`: Number of aromatic rings
- `fraction_csp3`: Fraction of sp³ hybridized carbons
- `qed`: Drug-likeness score

In addition to these properties, drug-likeness filters such as Lipinski's Rule of Five, Veber's Rules, and PAINS filter can also be applied.

## License

MIT