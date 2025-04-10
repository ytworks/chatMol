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