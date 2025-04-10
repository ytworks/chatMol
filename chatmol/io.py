"""
Module for handling input/output of molecular data
"""
import io
import logging
import csv
from typing import Any, Dict, List, Optional, Union, Callable

import pandas as pd

from .properties import (calculate_properties, calculate_all_properties, 
                       calculate_selected_properties, calculate_molecular_features,
                       MOLECULAR_FILTERS)

# Logger configuration
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def process_csv_data(csv_content: str, smiles_column: Optional[str] = None, 
                   properties: List[str] = None, filters: List[str] = None) -> Dict[str, Any]:
    """
    Calculate molecular properties and apply druglikeness filters from SMILES columns in CSV data and add them as new columns
    
    Args:
        csv_content: CSV data content to process (text format)
        smiles_column: Column name containing SMILES structures (if omitted, uses the rightmost column)
        properties: List of properties to calculate (molecular_weight, logp,
                   num_h_donors, num_h_acceptors, formula)
        filters: List of druglikeness filters to apply ('lipinski', 'veber', 'ghose', 
                 'egan', 'muegge', 'pains', 'all')
    
    Returns:
        Dict: Processing results
    """
    try:
        # Set default values
        if properties is None:
            properties = ["molecular_weight"]
            
        if filters is None:
            filters = []
            
        if not csv_content:
            return {
                "error": "No CSV data provided"
            }
        
        # Parse CSV data
        csv_data = io.StringIO(csv_content)
        df = pd.read_csv(csv_data)
        
        # Identify SMILES column
        if not smiles_column:
            smiles_column = df.columns[-1]  # Default is the rightmost column
            
        if smiles_column not in df.columns:
            return {
                "error": f"Specified SMILES column '{smiles_column}' not found in CSV data"
            }
        
        # Prepare DataFrame for results
        result_df = df.copy()
        
        # 一度のループで各SMILESに対して物性値とフィルターを同時に計算
        for idx, row in result_df.iterrows():
            smiles = row[smiles_column]
            
            # 統合関数を使って一度にすべての計算を実行
            features_result = calculate_molecular_features(smiles, properties, filters)
            
            # 物性値の結果をDataFrameに追加
            for prop_name, prop_value in features_result["properties"].items():
                column_name = prop_name
                # カラム名が既に存在する場合は名前を変更
                if prop_name in df.columns and prop_name != smiles_column:
                    column_name = f"{prop_name}_calculated"
                
                # 列が存在しない場合は追加
                if column_name not in result_df.columns:
                    result_df[column_name] = None
                
                # 値を設定
                result_df.at[idx, column_name] = prop_value
            
            # フィルター判定結果をDataFrameに追加
            for filter_name in filters:
                if filter_name in MOLECULAR_FILTERS:
                    column_name = f"{filter_name}_pass"
                    
                    # 列が存在しない場合は追加
                    if column_name not in result_df.columns:
                        result_df[column_name] = None
                    
                    # 値を設定（統合関数の結果から直接取得）
                    if column_name in features_result:
                        result_df.at[idx, column_name] = features_result[column_name]
            
            # 全フィルター通過フラグを追加（フィルターが指定されている場合）
            if filters and "all_filters_passed" in features_result:
                if "all_filters_passed" not in result_df.columns:
                    result_df["all_filters_passed"] = None
                result_df.at[idx, "all_filters_passed"] = features_result["all_filters_passed"]
        
        # Return results in CSV format
        output = io.StringIO()
        result_df.to_csv(output, index=False)
        csv_result = output.getvalue()
        
        # 適用されたフィルター名を取得
        applied_filters = [f for f in filters if f in MOLECULAR_FILTERS]
        if 'all' in filters:
            applied_filters = list(MOLECULAR_FILTERS.keys())
        
        return {
            "result": csv_result,
            "message": f"Molecular properties and filters applied. Rows processed: {len(result_df)}",
            "smiles_column": smiles_column,
            "properties_added": properties,
            "filters_applied": applied_filters
        }
        
    except Exception as e:
        logger.exception("Error occurred during processing")
        return {
            "error": f"An error occurred: {str(e)}"
        }


def read_smiles_from_csv(file_path: str, smiles_column: Optional[str] = None) -> List[str]:
    """
    Read a list of SMILES strings from a CSV file
    
    Args:
        file_path: Path to CSV file
        smiles_column: Column name containing SMILES structures (if omitted, uses the rightmost column)
    
    Returns:
        List[str]: List of SMILES strings
    """
    try:
        df = pd.read_csv(file_path)
        
        # Identify SMILES column
        if not smiles_column:
            smiles_column = df.columns[-1]  # Default is the rightmost column
            
        if smiles_column not in df.columns:
            logger.error(f"Specified SMILES column '{smiles_column}' not found in CSV file")
            return []
            
        return df[smiles_column].tolist()
        
    except Exception as e:
        logger.error(f"Error occurred while reading CSV file: {str(e)}")
        return []


def write_results_to_csv(file_path: str, data: List[Dict[str, Any]]) -> bool:
    """
    Write calculation results to a CSV file
    
    Args:
        file_path: Path to output CSV file
        data: List of data to write (list of dictionaries)
    
    Returns:
        bool: True if write was successful
    """
    try:
        if not data:
            logger.warning("No data to write")
            return False
            
        # Convert data to DataFrame
        df = pd.DataFrame(data)
        
        # Output to CSV
        df.to_csv(file_path, index=False)
        logger.info(f"Data saved to {file_path}")
        return True
        
    except Exception as e:
        logger.error(f"Error occurred while writing to CSV file: {str(e)}")
        return False