#!/usr/bin/env python
"""
MCP server for molecular properties calculation
"""
import csv
import io
import json
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

# Logger configuration
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Check for required module imports
try:
    import pandas as pd
    from mcp.server.fastmcp import FastMCP
    
    # Import chatMol library
    from chatmol.properties import calculate_molecular_features, get_available_properties
    from chatmol.properties import get_feature_descriptions
    from chatmol.io import add_properties_to_dataframe
    
    rdkit_available = True
except ImportError as e:
    print(f"Failed to import required modules: {str(e)}", file=sys.stderr)
    print("Please install required packages: pip install pandas mcp[server] rdkit", file=sys.stderr)
    # Implement minimal MCP server
    if 'mcp.server.fastmcp' in str(e):
        print("MCP module is not installed. Starting minimal server.", file=sys.stderr)
        try:
            # Minimal implementation using standard library only
            print("Starting minimal MCP server...", file=sys.stderr)
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
                        "message": "Required packages are not installed. Please run `pip install pandas mcp[server] rdkit`."
                    }
                    self.wfile.write(json.dumps(response).encode())
                
                def do_POST(self):
                    self.send_response(200)
                    self.send_header('Content-type', 'application/json')
                    self.end_headers()
                    response = {
                        "error": "Required packages are not installed. Please run `pip install pandas mcp[server] rdkit`."
                    }
                    self.wfile.write(json.dumps(response).encode())
            
            print(f"Starting minimal server on port {PORT}...", file=sys.stderr)
            with socketserver.TCPServer(("", PORT), MinimalHandler) as httpd:
                print(f"Server is running on port {PORT}", file=sys.stderr)
                httpd.serve_forever()
        except Exception as server_error:
            print(f"Failed to start minimal server: {str(server_error)}", file=sys.stderr)
            sys.exit(1)
    sys.exit(1)

# Initialize MCP server
mcp = FastMCP("Molecular Properties Calculator")


@mcp.tool()
def calculate_molecular_properties(
    input_data: str, 
    input_type: str = "smiles", 
    smiles_column: Optional[str] = None
) -> Dict[str, Any]:
    """
    Calculate molecular properties for SMILES strings or CSV data
    
    Args:
        input_data: Either a single SMILES string or CSV data content
        input_type: Type of input data - "smiles" for a single SMILES string or "csv" for CSV data
        smiles_column: Column name containing SMILES structures (for CSV input, if omitted, uses the rightmost column)
    
    Returns:
        Dict: Dictionary containing calculated molecular properties
    """
    try:
        # If input data is empty
        if not input_data:
            return {"error": "No input data provided"}
        
        # Processing single SMILES
        if input_type.lower() == "smiles":
            # Calculate and return directly for a single SMILES
            features = calculate_molecular_features(input_data)
            return features
                
        # Processing CSV format        
        elif input_type.lower() == "csv":
            import io
            import pandas as pd
            import os
            
            # Determine if input is file path or CSV data and process accordingly
            if os.path.exists(input_data) and input_data.lower().endswith('.csv'):
                # Process as file path
                try:
                    df = pd.read_csv(input_data)
                    logger.info(f"CSV file loaded successfully from path: {input_data}")
                except Exception as e:
                    return {"error": f"Failed to read CSV file from path {input_data}: {str(e)}"}
            else:
                # Process as CSV data string
                try:
                    # Convert potential string line breaks (\\n) to actual line breaks
                    formatted_input = input_data.replace('\\n', '\n')
                    csv_data = io.StringIO(formatted_input)
                    df = pd.read_csv(csv_data)
                    logger.info("CSV data parsed successfully from string input")
                except Exception as e:
                    return {"error": f"Failed to parse CSV data from string: {str(e)}"}
                
            # Identify SMILES column
            if not smiles_column:
                smiles_column = df.columns[-1]  # Default is rightmost column
                
            if smiles_column not in df.columns:
                return {
                    "error": f"Specified SMILES column '{smiles_column}' not found in CSV data. Available columns: {', '.join(df.columns)}"
                }
            
            # DataFrame to store results
            result_df = df.copy()
            
            # Process all SMILES at once
            smiles_list = result_df[smiles_column].tolist()
            feature_results = []
            
            # Calculate properties for each SMILES
            for smiles in smiles_list:
                if pd.isna(smiles):  # Check for missing values
                    feature_results.append({"error": "Invalid or missing SMILES"})
                    continue
                    
                try:
                    features = calculate_molecular_features(smiles)
                    feature_results.append(features)
                except Exception as e:
                    feature_results.append({"error": f"Error processing {smiles}: {str(e)}"})
            
            # Add properties to results
            add_properties_to_dataframe(result_df, feature_results)
            
            # Output in CSV format
            output = io.StringIO()
            result_df.to_csv(output, index=False)
            csv_result = output.getvalue()
            
            return {
                "result_format": "csv",
                "result": csv_result,
                "message": f"Processed {len(smiles_list)} compounds"
            }
        else:
            return {"error": f"Unsupported input_type: {input_type}. Use 'smiles' or 'csv'."}
            
    except Exception as e:
        logger.exception("Error occurred during property calculation")
        return {"error": f"An error occurred: {str(e)}"}


@mcp.tool()
def get_available_features() -> Dict[str, Any]:
    """
    Get a list of all available molecular features that can be calculated
    
    Returns:
        Dict: Dictionary containing lists of all available properties and filters
    """
    try:
        # Get list of available properties
        properties = get_available_properties()
        
        # Use new get_feature_descriptions
        feature_descriptions = get_feature_descriptions()
        
        # Separate properties and filters
        filters = [name for name, info in feature_descriptions.items() 
                  if info.get("is_filter", False)]
        
        # Get detailed information for properties
        property_info = {}
        for prop in properties:
            if prop in feature_descriptions:
                property_info[prop] = {
                    "name_ja": feature_descriptions[prop].get("ja", ""),
                    "name_en": feature_descriptions[prop].get("en", ""),
                    "description": feature_descriptions[prop].get("description", "")
                }
            else:
                property_info[prop] = {
                    "name_ja": prop,
                    "name_en": prop,
                    "description": ""
                }
        
        # Get detailed information for filters
        filter_info = {}
        for filter_name in filters:
            if filter_name in feature_descriptions:
                filter_info[filter_name] = {
                    "description": feature_descriptions[filter_name].get("description", ""),
                    "result_key": feature_descriptions[filter_name].get("result_key", "")
                }
        
        return {
            "properties": properties,
            "property_count": len(properties),
            "property_details": property_info,
            "filters": filters,
            "filter_count": len(filters),
            "filter_details": filter_info,
            "message": "Available molecular features that can be calculated"
        }
        
    except Exception as e:
        logger.exception("Error occurred while retrieving available features")
        return {"error": f"An error occurred: {str(e)}"}


if __name__ == "__main__":
    try:
        # Check if RDKit is available
        mol_features = calculate_molecular_features("C")
        if not isinstance(mol_features, dict):
            print("RDKit does not appear to be properly installed", file=sys.stderr)
            sys.exit(1)
            
        print("Starting MCP server for molecular properties calculation...", file=sys.stderr)
        mcp.run()
    except Exception as e:
        print(f"Server startup error: {str(e)}", file=sys.stderr)
        sys.exit(1)