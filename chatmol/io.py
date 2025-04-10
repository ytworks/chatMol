"""
分子データの入出力を処理するモジュール
"""
import io
import logging
import csv
from typing import Any, Dict, List, Optional, Union

import pandas as pd

from .properties import calculate_properties

# ロガーの設定
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def process_csv_data(csv_content: str, smiles_column: Optional[str] = None, 
                   properties: List[str] = None) -> Dict[str, Any]:
    """
    CSVデータ内のSMILES列から分子特性を計算し、新しい列として追加します
    
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


def read_smiles_from_csv(file_path: str, smiles_column: Optional[str] = None) -> List[str]:
    """
    CSVファイルからSMILES文字列のリストを読み込む
    
    Args:
        file_path: CSVファイルのパス
        smiles_column: SMILES構造を含む列名（指定しない場合は一番右の列を使用）
    
    Returns:
        List[str]: SMILES文字列のリスト
    """
    try:
        df = pd.read_csv(file_path)
        
        # SMILES列を特定
        if not smiles_column:
            smiles_column = df.columns[-1]  # デフォルトは一番右の列
            
        if smiles_column not in df.columns:
            logger.error(f"指定されたSMILES列 '{smiles_column}' がCSVファイルに見つかりません")
            return []
            
        return df[smiles_column].tolist()
        
    except Exception as e:
        logger.error(f"CSVファイルの読み込み中にエラーが発生しました: {str(e)}")
        return []


def write_results_to_csv(file_path: str, data: List[Dict[str, Any]]) -> bool:
    """
    計算結果をCSVファイルに書き込む
    
    Args:
        file_path: 出力先CSVファイルのパス
        data: 書き込むデータのリスト（辞書のリスト形式）
    
    Returns:
        bool: 書き込みが成功した場合はTrue
    """
    try:
        if not data:
            logger.warning("書き込むデータがありません")
            return False
            
        # データをDataFrameに変換
        df = pd.DataFrame(data)
        
        # CSVに出力
        df.to_csv(file_path, index=False)
        logger.info(f"データを {file_path} に保存しました")
        return True
        
    except Exception as e:
        logger.error(f"CSVファイルへの書き込み中にエラーが発生しました: {str(e)}")
        return False