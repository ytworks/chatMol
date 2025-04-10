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
    Add molecular property calculation results in flat format to a DataFrame
    
    Args:
        df: DataFrame to add properties to
        feature_results: List of molecular property calculation results in flat format
        
    Returns:
        None: DataFrame is updated by reference
    """
    # Get all keys
    all_keys = set()
    for result in feature_results:
        all_keys.update(result.keys())
    
    # Exclude specific keys (smiles, error, mol, etc.)
    exclude_keys = {"smiles", "error", "mol", "pains_alerts"}
    properties = [key for key in all_keys if key not in exclude_keys]
    
    # Add each property to the DataFrame
    for prop_name in properties:
        # Change column name if it conflicts with an existing one
        column_name = prop_name
        if prop_name in df.columns:
            column_name = f"{prop_name}_calculated"
            
        # Get values for each row
        values = []
        for result in feature_results:
            values.append(result.get(prop_name))
            
        # Add column to DataFrame
        df[column_name] = values