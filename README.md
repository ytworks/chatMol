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