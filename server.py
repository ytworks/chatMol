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
    format="%(asctime)s - %(name)s - %(levellevel)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Check for required module imports
try:
    import pandas as pd
    from mcp.server.fastmcp import FastMCP
    
    # Import chatMol library
    from chatmol.properties import calculate_molecular_weight, calculate_properties
    from chatmol.io import process_csv_data
    
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
mcp = FastMCP("Molecular Weight Calculator")


@mcp.tool()
def add_molecular_weight(csv_content: str, smiles_column: Optional[str] = None, properties: List[str] = None) -> Dict[str, Any]:
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
    # Use process_csv_data function from chatMol library
    return process_csv_data(csv_content, smiles_column, properties)


if __name__ == "__main__":
    try:
        # Check if RDKit is available by using the chatMol library function
        result = calculate_molecular_weight("C")
        if not isinstance(result, float):
            print("RDKit does not appear to be properly installed", file=sys.stderr)
            sys.exit(1)
            
        print("Starting MCP server...", file=sys.stderr)
        mcp.run()
    except Exception as e:
        print(f"Server startup error: {str(e)}", file=sys.stderr)
        sys.exit(1)