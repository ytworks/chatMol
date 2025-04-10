#!/usr/bin/env python
"""
分子量計算のためのMCPサーバー
"""
import csv
import io
import json
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

# ロガーの設定
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# 必須モジュールのインポートチェック
try:
    import pandas as pd
    from mcp.server.fastmcp import FastMCP
    rdkit_available = True
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
    except ImportError:
        print("RDKitモジュールをインポートできません。システムにインストールされているか確認してください。", file=sys.stderr)
        rdkit_available = False
except ImportError as e:
    print(f"必須モジュールのインポートに失敗しました: {str(e)}", file=sys.stderr)
    print("必須パッケージをインストールしてください: pip install pandas mcp[server] rdkit", file=sys.stderr)
    # 最小限のMCPサーバーを実装
    if 'mcp.server.fastmcp' in str(e):
        print("MCPモジュールがインストールされていません。最小限のサーバーを起動します。", file=sys.stderr)
        try:
            # 標準ライブラリのみを使用した最小限実装
            print("最小限のMCPサーバーを起動中...", file=sys.stderr)
            import http.server
            import socketserver
            
            PORT = 8080
            
            class MinimalHandler(http.server.SimpleHTTPRequestHandler):
                def do_GET(self):
                    self.send_response(200)
                    self.send_header('Content-type', 'application/json')
                    self.end_headers()
                    response = {
                        "status": "error",
                        "message": "必須パッケージがインストールされていません。`pip install pandas mcp[server] rdkit`を実行してください。"
                    }
                    self.wfile.write(json.dumps(response).encode())
                
                def do_POST(self):
                    self.send_response(200)
                    self.send_header('Content-type', 'application/json')
                    self.end_headers()
                    response = {
                        "error": "必須パッケージがインストールされていません。`pip install pandas mcp[server] rdkit`を実行してください。"
                    }
                    self.wfile.write(json.dumps(response).encode())
            
            print(f"最小限のサーバーをポート {PORT} で起動中...", file=sys.stderr)
            with socketserver.TCPServer(("", PORT), MinimalHandler) as httpd:
                print(f"サーバーがポート {PORT} で実行中", file=sys.stderr)
                httpd.serve_forever()
        except Exception as server_error:
            print(f"最小限サーバーの起動に失敗しました: {str(server_error)}", file=sys.stderr)
            sys.exit(1)
    sys.exit(1)

# MCPサーバーの初期化
mcp = FastMCP("Molecular Weight Calculator")


def calculate_molecular_weight(smiles: str) -> float:
    """
    SMILES文字列から分子量を計算する

    Args:
        smiles: SMILES表記の分子構造

    Returns:
        float: 分子量。SMILES文字列が無効の場合はNaN
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return float('nan')
        return Descriptors.MolWt(mol)
    except Exception:
        return float('nan')


def calculate_properties(smiles: str) -> Dict[str, Union[float, str, None]]:
    """
    SMILES文字列から分子の基本的な特性を取得

    Args:
        smiles: SMILES表記の分子構造

    Returns:
        Dict: 分子の基本特性を含む辞書
    """
    props = {
        "molecular_weight": float('nan'),
        "logp": float('nan'),
        "num_h_donors": float('nan'),
        "num_h_acceptors": float('nan'),
        "formula": None
    }
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return props
            
        props["molecular_weight"] = Descriptors.MolWt(mol)
        props["logp"] = Descriptors.MolLogP(mol)
        props["num_h_donors"] = Descriptors.NumHDonors(mol)
        props["num_h_acceptors"] = Descriptors.NumHAcceptors(mol)
        props["formula"] = Chem.rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        pass
        
    return props


@mcp.tool()
def add_molecular_weight(csv_content: str, smiles_column: Optional[str] = None, properties: List[str] = None) -> Dict[str, Any]:
    """
    CSVデータ内のSMILES列から分子量などの特性を計算し、新しい列として追加します
    
    Args:
        csv_content: 処理するCSVデータの内容（テキスト形式）
        smiles_column: SMILES構造を含む列名（指定しない場合は一番右の列を使用）
        properties: 計算する特性のリスト（分子量:molecular_weight、脂溶性:logp、
                   水素結合ドナー数:num_h_donors、水素結合アクセプター数:num_h_acceptors、分子式:formula）
    
    Returns:
        Dict: 処理結果
    """
    try:
        # デフォルト値の設定
        if properties is None:
            properties = ["molecular_weight"]
            
        if not csv_content:
            return {
                "error": "CSVデータが提供されていません"
            }
        
        # CSVデータをパース
        csv_data = io.StringIO(csv_content)
        df = pd.read_csv(csv_data)
        
        # SMILES列を特定
        if not smiles_column:
            smiles_column = df.columns[-1]  # デフォルトは一番右の列
            
        if smiles_column not in df.columns:
            return {
                "error": f"指定されたSMILES列 '{smiles_column}' がCSVデータに見つかりません"
            }
        
        # 結果用のDataFrameを準備
        result_df = df.copy()
        
        # 指定されたプロパティを計算して追加
        for smiles_idx, smiles in enumerate(result_df[smiles_column]):
            props = calculate_properties(smiles)
            
            # 要求されたプロパティのみ追加
            for prop_name in properties:
                if prop_name in props:
                    if result_df.shape[0] > smiles_idx:  # インデックスが範囲内か確認
                        column_name = prop_name
                        # 既に同名の列が存在する場合は名前を変更
                        if prop_name in result_df.columns and prop_name != smiles_column:
                            column_name = f"{prop_name}_calculated"
                        
                        result_df.loc[smiles_idx, column_name] = props[prop_name]
        
        # 結果をCSV形式で返す
        output = io.StringIO()
        result_df.to_csv(output, index=False)
        csv_result = output.getvalue()
        
        return {
            "result": csv_result,
            "message": f"分子特性の計算が完了しました。処理された行数: {len(result_df)}",
            "smiles_column": smiles_column,
            "properties_added": properties
        }
        
    except Exception as e:
        logger.exception("処理中にエラーが発生しました")
        return {
            "error": f"エラーが発生しました: {str(e)}"
        }


if __name__ == "__main__":
    if rdkit_available:
        print("MCPサーバーを起動します...", file=sys.stderr)
        try:
            mcp.run()
        except Exception as e:
            print(f"サーバー起動エラー: {str(e)}", file=sys.stderr)
            sys.exit(1)
    else:
        print("RDKitが利用できないため、サーバーを起動できません", file=sys.stderr)
        sys.exit(1)