"""
Module for handling input/output of molecular data
"""
import io
import logging
import csv
from typing import Any, Dict, List, Optional, Union

import pandas as pd

from .properties import calculate_properties

# Logger configuration
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def process_csv_data(csv_content: str, smiles_column: Optional[str] = None, 
                   properties: List[str] = None) -> Dict[str, Any]:
    """
    Calculate molecular properties from SMILES columns in CSV data and add them as new columns
    
    Args:
        csv_content: CSV data content to process (text format)
        smiles_column: Column name containing SMILES structures (if omitted, uses the rightmost column)
        properties: List of properties to calculate (molecular_weight, logp,
                   num_h_donors, num_h_acceptors, formula)
    
    Returns:
        Dict: Processing results
    """
    try:
        # Set default values
        if properties is None:
            properties = ["molecular_weight"]
            
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
        
        # Initialize the property columns
        for prop in properties:
            result_df[prop] = None
        
        # Calculate and add specified properties
        for idx, row in result_df.iterrows():
            smiles = row[smiles_column]
            props = calculate_properties(smiles)
            
            # Add only requested properties
            for prop_name in properties:
                if prop_name in props:
                    column_name = prop_name
                    # Rename column if it already exists
                    if prop_name in df.columns and prop_name != smiles_column:
                        column_name = f"{prop_name}_calculated"
                    
                    result_df.at[idx, column_name] = props[prop_name]
        
        # Return results in CSV format
        output = io.StringIO()
        result_df.to_csv(output, index=False)
        csv_result = output.getvalue()
        
        return {
            "result": csv_result,
            "message": f"Molecular property calculation completed. Rows processed: {len(result_df)}",
            "smiles_column": smiles_column,
            "properties_added": properties
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