#!/usr/bin/env python
"""
MCP server for molecular weight calculation
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
    from chatmol.properties import (
        calculate_molecular_weight, calculate_properties, 
        calculate_all_properties, calculate_selected_properties,
        get_available_properties, get_properties_table,
        get_property_descriptions, get_property_categories,
        check_lipinski_rule, check_veber_rules, check_ghose_filter,
        check_egan_filter, check_muegge_filter, check_pains_filter,
        get_druglikeness_filters_summary, check_all_druglikeness_filters
    )
    from chatmol.io import process_csv_data, DRUGLIKENESS_FILTERS
    
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
def add_molecular_weight(csv_content: str, smiles_column: Optional[str] = None, properties: List[str] = None, filters: List[str] = None) -> Dict[str, Any]:
    """
    Calculate molecular properties and apply druglikeness filters from SMILES columns in CSV data
    
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
    # Use process_csv_data function from chatMol library
    return process_csv_data(csv_content, smiles_column, properties, filters)


@mcp.tool()
def calculate_molecular_properties(smiles: str, properties: List[str] = None) -> Dict[str, Any]:
    """
    Calculate molecular properties for a single SMILES string
    
    Args:
        smiles: Molecular structure in SMILES notation
        properties: List of specific properties to calculate. If None, returns basic properties only.
                   If ['all'], returns all available properties.
    
    Returns:
        Dict: Dictionary containing calculated molecular properties
    """
    try:
        if not smiles:
            return {"error": "No SMILES string provided"}
        
        # 基本的な物性値を計算（デフォルト）
        if properties is None:
            return calculate_properties(smiles)
        
        # 全ての物性値を計算
        if len(properties) == 1 and properties[0].lower() == 'all':
            return calculate_all_properties(smiles)
        
        # 指定された物性値のみを計算
        return calculate_selected_properties(smiles, properties)
        
    except Exception as e:
        logger.exception("Error occurred during property calculation")
        return {"error": f"An error occurred: {str(e)}"}


@mcp.tool()
def get_property_list() -> Dict[str, Any]:
    """
    Get a list of all available molecular properties that can be calculated
    
    Returns:
        Dict: Dictionary containing a list of all available property names
    """
    try:
        properties = get_available_properties()
        return {
            "available_properties": properties,
            "count": len(properties)
        }
    except Exception as e:
        logger.exception("Error occurred while retrieving property list")
        return {"error": f"An error occurred: {str(e)}"}


@mcp.tool()
def get_properties_info(format: str = "markdown", language: str = "ja", include_categories: bool = False) -> Dict[str, Any]:
    """
    Get formatted information about all available molecular properties
    
    Args:
        format: Output format ('markdown', 'csv', 'json', 'text')
        language: Language for property names ('ja' for Japanese, 'en' for English)
        include_categories: Whether to include property categories in the result
        
    Returns:
        Dict: Dictionary containing formatted property information
    """
    try:
        # 物性値の説明テーブルを取得
        properties_table = get_properties_table(format, language)
        
        result = {
            "properties_table": properties_table,
            "format": format,
            "language": language
        }
        
        # カテゴリ情報も要求された場合
        if include_categories:
            categories = get_property_categories()
            result["categories"] = categories
            
            # カテゴリ毎の物性値名を日本語/英語で提供
            if format == "json":
                property_info = get_property_descriptions()
                categorized_properties = {}
                for category, props in categories.items():
                    categorized_properties[category] = {
                        "properties": props,
                        "names": {
                            prop: {
                                "ja": property_info[prop]["ja"], 
                                "en": property_info[prop]["en"]
                            } for prop in props if prop in property_info
                        }
                    }
                result["categorized_properties"] = categorized_properties
        
        return result
        
    except Exception as e:
        logger.exception("Error occurred while retrieving property information")
        return {"error": f"An error occurred: {str(e)}"}


@mcp.tool()
def check_drug_likeness(smiles: str, rules: List[str] = None) -> Dict[str, Any]:
    """
    Check if a molecule complies with various drug-likeness rules
    
    Args:
        smiles: Molecular structure in SMILES notation
        rules: List of rule sets to check. Options: 'lipinski', 'veber', 'ghose', 
               'egan', 'muegge', 'pains', 'all'. Default is ['lipinski'] if not specified.
    
    Returns:
        Dict: Results for each requested rule check
    """
    try:
        if not smiles:
            return {"error": "No SMILES string provided"}
        
        # デフォルトでLipinski's Rule of Fiveをチェック
        if rules is None:
            rules = ["lipinski"]
            
        # 全てのルールをチェックするリクエスト
        if "all" in [r.lower() for r in rules]:
            # 全てのフィルターを実行
            full_results = check_all_druglikeness_filters(smiles)
            # 簡易概要も追加
            summary = get_druglikeness_filters_summary(smiles)
            
            # 分子のSMILESと基本物性を追加
            full_results["smiles"] = smiles
            full_results["basic_properties"] = calculate_properties(smiles)
            full_results["summary"] = summary
            
            return full_results
            
        result = {}
            
        # リクエストに応じて各ルールをチェック
        for rule in rules:
            rule_lower = rule.lower()
            if rule_lower == "lipinski":
                result["lipinski"] = check_lipinski_rule(smiles)
            elif rule_lower == "veber":
                result["veber"] = check_veber_rules(smiles)
            elif rule_lower == "ghose":
                result["ghose"] = check_ghose_filter(smiles)
            elif rule_lower == "egan":
                result["egan"] = check_egan_filter(smiles)
            elif rule_lower == "muegge":
                result["muegge"] = check_muegge_filter(smiles)
            elif rule_lower == "pains":
                result["pains"] = check_pains_filter(smiles)
                
        # 分子のSMILESと基本物性も返す
        result["smiles"] = smiles
        result["basic_properties"] = calculate_properties(smiles)
            
        return result
        
    except Exception as e:
        logger.exception("Error occurred during drug-likeness check")
        return {"error": f"An error occurred: {str(e)}"}


if __name__ == "__main__":
    try:
        # Check if RDKit is available by using the chatMol library function
        result = calculate_molecular_weight("C")
        if not isinstance(result, float):
            print("RDKit does not appear to be properly installed", file=sys.stderr)
            sys.exit(1)
            
        print("Starting MCP server for molecular properties calculation...", file=sys.stderr)
        mcp.run()
    except Exception as e:
        print(f"Server startup error: {str(e)}", file=sys.stderr)
        sys.exit(1)