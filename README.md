# Molecular Weight Calculator MCP

An Model Context Protocol (MCP) server for calculating molecular weights of compounds in CSV files from SMILES notation

## Overview

This MCP server calculates molecular weights and other molecular properties (lipophilicity, hydrogen bond donors, hydrogen bond acceptors, molecular formula, etc.) from SMILES columns in CSV files and adds the results as new columns. When used with Claude Desktop, it enables simple analysis of chemical data.

## Requirements

- Python 3.10 or higher
- RDKit 2024.9.6
- pandas 2.2.3
- mcp (Model Context Protocol) package v1.2.0 or higher
- Claude Desktop
- uv (Python package manager)

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
    "Molecular Weight Calculator": {
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

Note:
- If the `uv` command is not in your environment path, use an absolute path (e.g.: `/path/to/uv`).
- Replace `/path/to/chatMol/server.py` with the absolute path to this script.
- Always use absolute paths, not relative paths.

### 2. Available Tools

This MCP server provides the following tool:

#### add_molecular_weight

- Description: Calculates molecular properties from SMILES columns in CSV data and adds the results as new columns
- Input parameters:
  - `csv_content`: CSV data content to process (required)
  - `smiles_column`: Column name containing SMILES structures (if omitted, uses the rightmost column)
  - `properties`: List of properties to calculate (if omitted, calculates molecular weight only)
    - `molecular_weight`: Molecular weight
    - `logp`: Lipophilicity
    - `num_h_donors`: Number of hydrogen bond donors
    - `num_h_acceptors`: Number of hydrogen bond acceptors
    - `formula`: Molecular formula

### Examples

In Claude Desktop, you can calculate molecular weights from SMILES columns in CSV data as follows:

```
Please calculate the molecular weights for this CSV file:

ID,Name,SMILES
1,Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O
2,Paracetamol,CC(=O)NC1=CC=C(C=C1)O
3,Ibuprofen,CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
```

To calculate multiple properties simultaneously:

```
Please calculate molecular weight, logP, and molecular formula from the SMILES column in this CSV:

ID,Name,SMILES
1,Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O
2,Paracetamol,CC(=O)NC1=CC=C(C=C1)O
3,Ibuprofen,CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
```

You can also specify the column name containing SMILES:

```
Please calculate molecular weights from the smiles_col column in this CSV data:

ID,Name,smiles_col
1,Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O
2,Paracetamol,CC(=O)NC1=CC=C(C=C1)O
3,Ibuprofen,CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
```

## License

MIT