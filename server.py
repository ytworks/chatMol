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
    
    # chatMolライブラリをインポート
    from chatmol.properties import calculate_molecular_weight, calculate_properties
    from chatmol.io import process_csv_data
    
    rdkit_available = True
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
    # chatMolライブラリのprocess_csv_data関数を使用
    return process_csv_data(csv_content, smiles_column, properties)


if __name__ == "__main__":
    try:
        # RDKitが利用可能かどうかを確認するため、chatMolライブラリの関数を使用
        result = calculate_molecular_weight("C")
        if not isinstance(result, float):
            print("RDKitが正しくインストールされていないようです", file=sys.stderr)
            sys.exit(1)
            
        print("MCPサーバーを起動します...", file=sys.stderr)
        mcp.run()
    except Exception as e:
        print(f"サーバー起動エラー: {str(e)}", file=sys.stderr)
        sys.exit(1)