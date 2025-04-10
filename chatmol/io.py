"""
Module for handling input/output of molecular data
"""
import io
import logging
import csv
from typing import Any, Dict, List, Optional, Union, Callable

import pandas as pd

from .properties import (calculate_properties, check_lipinski_rule, check_veber_rules,
                         check_ghose_filter, check_egan_filter, check_muegge_filter,
                         check_pains_filter, get_druglikeness_filters_summary)

# Logger configuration
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# ドラッグライクネスフィルター関数と列名のマッピング
DRUGLIKENESS_FILTERS = {
    "lipinski": {
        "function": check_lipinski_rule,
        "result_key": "all_rules_passed",
        "column_name": "lipinski_pass"
    },
    "veber": {
        "function": check_veber_rules,
        "result_key": "all_rules_passed",
        "column_name": "veber_pass"
    },
    "ghose": {
        "function": check_ghose_filter,
        "result_key": "all_rules_passed",
        "column_name": "ghose_pass"
    },
    "egan": {
        "function": check_egan_filter,
        "result_key": "all_rules_passed",
        "column_name": "egan_pass"
    },
    "muegge": {
        "function": check_muegge_filter,
        "result_key": "all_rules_passed",
        "column_name": "muegge_pass"
    },
    "pains": {
        "function": check_pains_filter,
        "result_key": "pains_free",
        "column_name": "pains_free"
    }
}

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
        
        # 'all'フィルターが指定された場合、すべてのフィルターを適用
        if 'all' in filters:
            filters = list(DRUGLIKENESS_FILTERS.keys())
        
        # 適用するフィルターの情報を整理
        applied_filters = []
        filter_configs = []
        for filter_name in filters:
            if filter_name in DRUGLIKENESS_FILTERS:
                filter_info = DRUGLIKENESS_FILTERS[filter_name]
                column_name = filter_info["column_name"]
                result_df[column_name] = None
                applied_filters.append(filter_name)
                filter_configs.append({
                    "name": filter_name,
                    "function": filter_info["function"],
                    "result_key": filter_info["result_key"],
                    "column_name": column_name
                })
        
        # 計算するプロパティのカラムを初期化
        property_columns = {}
        for prop in properties:
            column_name = prop
            # カラム名が既に存在する場合は名前を変更
            if prop in df.columns and prop != smiles_column:
                column_name = f"{prop}_calculated"
            result_df[column_name] = None
            property_columns[prop] = column_name
        
        # 一度のループで各SMILESに対してプロパティとフィルターを計算
        for idx, row in result_df.iterrows():
            smiles = row[smiles_column]
            
            # プロパティを計算（必要な場合のみ）
            if properties:
                props = calculate_properties(smiles)
                for prop_name, column_name in property_columns.items():
                    if prop_name in props:
                        result_df.at[idx, column_name] = props[prop_name]
            
            # フィルターを適用（必要な場合のみ）
            for config in filter_configs:
                filter_func = config["function"]
                result_key = config["result_key"]
                column_name = config["column_name"]
                
                # フィルター関数を実行し、結果を取得
                filter_result = filter_func(smiles)
                
                # ブール値結果を取り出して、DataFrame列に追加
                result_df.at[idx, column_name] = filter_result[result_key]
        
        # Return results in CSV format
        output = io.StringIO()
        result_df.to_csv(output, index=False)
        csv_result = output.getvalue()
        
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