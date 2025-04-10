"""
Module for handling input/output of molecular data
"""
import io
import logging
from typing import Any, Dict, List

import pandas as pd

from .properties import calculate_molecular_features

# Logger configuration
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def convert_properties_to_markdown(smiles: str, features: Dict[str, Any]) -> str:
    """
    単一分子のプロパティをマークダウン形式に変換する
    
    Args:
        smiles: 分子のSMILES表記
        features: フラット形式の分子特性辞書
        
    Returns:
        str: マークダウン形式の文字列
    """
    md_lines = []
    
    # 分子情報のヘッダー
    md_lines.append(f"## 分子構造情報: {smiles}")
    md_lines.append("")
    
    # 一般的な物性値
    md_lines.append("### 分子特性")
    md_lines.append("")
    md_lines.append("| プロパティ | 値 |")
    md_lines.append("|-----------|-----|")
    
    # フィルター関連のキーを除外
    for prop, value in sorted(features.items()):
        # フィルター関連、特殊キーは除外
        if (not prop.endswith("_pass") and not prop.endswith("_ok") and 
            prop != "smiles" and prop != "pains_free" and 
            prop != "all_filters_passed" and prop != "pains_alerts" and
            prop != "pains_num_alerts" and prop != "error" and prop != "mol"):
            
            # 値の整形（数値の場合は丸める）
            if isinstance(value, float):
                formatted_value = f"{value:.3f}"
            else:
                formatted_value = str(value) if value is not None else "N/A"
                
            md_lines.append(f"| {prop} | {formatted_value} |")
    
    md_lines.append("")
    
    # フィルター結果
    md_lines.append("### ドラッグライクネスフィルター")
    md_lines.append("")
    
    # Lipinski
    if "lipinski_pass" in features:
        md_lines.append("#### Lipinski")
        md_lines.append("")
        md_lines.append("| 条件 | 結果 |")
        md_lines.append("|------|------|")
        for key in ["lipinski_molecular_weight_ok", "lipinski_logp_ok", "lipinski_h_donors_ok", "lipinski_h_acceptors_ok"]:
            if key in features:
                result_str = "✓" if features.get(key) else "✗"
                md_lines.append(f"| {key.replace('lipinski_', '').replace('_ok', '')} | {result_str} |")
        md_lines.append(f"| **全ルール** | {'✓' if features.get('lipinski_pass') else '✗'} |")
        md_lines.append("")
    
    # 他のフィルター
    filters = ["veber", "ghose", "egan", "muegge"]
    for filter_name in filters:
        pass_key = f"{filter_name}_pass"
        if pass_key in features:
            md_lines.append(f"#### {filter_name.capitalize()}")
            md_lines.append("")
            md_lines.append("| 条件 | 結果 |")
            md_lines.append("|------|------|")
            
            # フィルター条件を探す
            filter_keys = [k for k in features.keys() if k.startswith(f"{filter_name}_") and k != pass_key]
            for key in filter_keys:
                result_str = "✓" if features.get(key) else "✗"
                md_lines.append(f"| {key.replace(f'{filter_name}_', '').replace('_ok', '')} | {result_str} |")
                
            md_lines.append(f"| **全ルール** | {'✓' if features.get(pass_key) else '✗'} |")
            md_lines.append("")
    
    # PAINSフィルター
    if "pains_free" in features:
        md_lines.append("#### PAINS")
        md_lines.append("")
        md_lines.append("| 条件 | 結果 |")
        md_lines.append("|------|------|")
        md_lines.append(f"| PAINS-free | {'✓' if features.get('pains_free') else '✗'} |")
        
        if not features.get("pains_free") and "pains_alerts" in features and features["pains_alerts"]:
            md_lines.append("")
            md_lines.append("検出されたPAINSパターン:")
            for alert in features["pains_alerts"]:
                md_lines.append(f"- {alert.get('description', 'Unknown')}")
        md_lines.append("")
    
    return "\n".join(md_lines)


def convert_dataframe_to_markdown(df: pd.DataFrame) -> str:
    """
    DataFrame全体をマークダウンテーブル形式に変換する
    
    Args:
        df: 変換するDataFrame
        
    Returns:
        str: マークダウン形式のテーブル
    """
    # pandas DataFrameのto_markdownメソッドを使用
    try:
        return df.to_markdown(index=False)
    except AttributeError:
        # to_markdownが利用できない場合、手動で変換
        header = "| " + " | ".join(df.columns) + " |"
        separator = "| " + " | ".join(["---"] * len(df.columns)) + " |"
        
        rows = []
        for _, row in df.iterrows():
            formatted_row = []
            for value in row:
                if isinstance(value, float):
                    formatted_row.append(f"{value:.3f}")
                else:
                    formatted_row.append(str(value) if value is not None else "N/A")
            rows.append("| " + " | ".join(formatted_row) + " |")
        
        return "\n".join([header, separator] + rows)


def add_properties_to_dataframe(df: pd.DataFrame, feature_results: List[Dict[str, Any]]) -> None:
    """
    フラット形式の分子特性計算結果をDataFrameに追加する
    
    Args:
        df: 特性を追加するDataFrame
        feature_results: フラット形式の分子特性計算結果のリスト
        
    Returns:
        None: DataFrameは参照で更新
    """
    # すべてのキーを取得
    all_keys = set()
    for result in feature_results:
        all_keys.update(result.keys())
    
    # 特定のキーを除外（smiles, error, mol など）
    exclude_keys = {"smiles", "error", "mol", "pains_alerts"}
    properties = [key for key in all_keys if key not in exclude_keys]
    
    # 各プロパティをDataFrameに追加
    for prop_name in properties:
        # カラム名が既存のものとぶつかる場合は名前を変更
        column_name = prop_name
        if prop_name in df.columns:
            column_name = f"{prop_name}_calculated"
            
        # 各行の値を取得
        values = []
        for result in feature_results:
            values.append(result.get(prop_name))
            
        # DataFrameにカラムを追加
        df[column_name] = values