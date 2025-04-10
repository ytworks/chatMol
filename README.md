# Molecular Weight Calculator MCP

分子構造のSMILESからCSVファイル内の化合物の分子量を計算するためのModel Context Protocol (MCP) サーバー

## 概要

このMCPサーバーは、CSVファイル内のSMILES列から分子量や他の分子特性（脂溶性、水素結合ドナー数、水素結合アクセプター数、分子式など）を計算し、結果を新しい列として追加します。Claude Desktopと連携して使用することで、化学データの簡単な解析が可能になります。

## 環境要件

- Python 3.10以上
- RDKit 2024.9.6
- pandas 2.2.3
- mcp (Model Context Protocol) パッケージ v1.2.0以上
- Claude Desktop
- uv (Pythonパッケージマネージャー)

## インストール方法

### 1. 依存関係のインストール

```bash
# uvを使用する場合
uv pip install "rdkit==2024.9.6" "pandas==2.2.3" "mcp[cli,server]>=1.2.0"

# pipを使用する場合
pip install "rdkit==2024.9.6" "pandas==2.2.3" "mcp[cli,server]>=1.2.0"
```

## 使用方法

### 1. Claude Desktopでの使用

#### 設定ファイルの編集

Claude Desktopの設定ファイルを編集して、このMCPサーバーを追加します。設定ファイルのパスは以下の通りです：

- macOS: `~/Library/Application Support/Claude/claude_desktop_config.json`
- Windows: `%APPDATA%\Claude\claude_desktop_config.json`

以下のJSON設定を追加します（既存の`mcpServers`オブジェクト内に追加）：

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

注意:
- `uv` コマンドが環境変数のパスに含まれていない場合は、絶対パス（例：`/path/to/uv`）を使用してください。
- `/path/to/chatMol/server.py` の部分は、このスクリプトの絶対パスを指定してください。
- 相対パスは使用せず、必ず絶対パスを使用してください。

### 2. 使用可能なツール

このMCPサーバーは以下のツールを提供します：

#### add_molecular_weight

- 説明：CSVデータ内のSMILES列から分子量などの特性を計算し、結果を新しい列として追加します
- 入力パラメータ：
  - `csv_content`：処理するCSVデータの内容（必須）
  - `smiles_column`：SMILES構造を含む列名（省略時は一番右の列を使用）
  - `properties`：計算する特性のリスト（省略時は分子量のみ）
    - `molecular_weight`：分子量
    - `logp`：脂溶性
    - `num_h_donors`：水素結合ドナー数
    - `num_h_acceptors`：水素結合アクセプター数
    - `formula`：分子式

### 使用例

Claude Desktopで以下のように入力することで、CSVデータ内のSMILES列から分子量を計算できます：

```
このCSVファイルの分子量を計算してください：

ID,Name,SMILES
1,Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O
2,Paracetamol,CC(=O)NC1=CC=C(C=C1)O
3,Ibuprofen,CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
```

複数の特性を同時に計算する場合：

```
このCSVのSMILES列から分子量、脂溶性、分子式を計算してください：

ID,Name,SMILES
1,Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O
2,Paracetamol,CC(=O)NC1=CC=C(C=C1)O
3,Ibuprofen,CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
```

また、SMILES列の列名を指定することもできます：

```
このCSVデータのsmiles_col列から分子量を計算してください：

ID,Name,smiles_col
1,Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O
2,Paracetamol,CC(=O)NC1=CC=C(C=C1)O
3,Ibuprofen,CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
```

## ライセンス

MIT