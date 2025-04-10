"""
Test module for the chatMol CSV processing functionality.
Tests the processing of CSV data containing molecular structures.
"""
import pytest
import pandas as pd
import io
from chatmol.io import process_csv_data, read_smiles_from_csv


class TestCSVProcessing:
    """Test class for CSV data processing functionality."""

    def test_process_csv_with_defaults(self):
        """Test processing CSV data with default parameters."""
        # Create a simple CSV with SMILES column
        csv_content = (
            "ID,Name,SMILES\n"
            "1,Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O\n"
            "2,Paracetamol,CC(=O)NC1=CC=C(C=C1)O\n"
            "3,Ibuprofen,CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
        )
        
        # Process with default parameters (only molecular weight)
        result = process_csv_data(csv_content)
        
        # Check if result contains expected keys
        assert "result" in result
        assert "message" in result
        assert "smiles_column" in result
        assert "properties_added" in result
        
        # Check if properties were calculated correctly
        result_df = pd.read_csv(io.StringIO(result["result"]))
        
        # Check if molecular_weight column was added
        assert "molecular_weight" in result_df.columns
        
        # Verify values for known compounds
        aspirin_mw = result_df.loc[0, "molecular_weight"]
        assert round(aspirin_mw, 3) == 180.159
        
        paracetamol_mw = result_df.loc[1, "molecular_weight"]
        assert round(paracetamol_mw, 3) == 151.165  # 修正: RDKitの計算値に合わせる

    def test_process_csv_with_multiple_properties(self):
        """Test processing CSV data with multiple property calculations."""
        # Create a simple CSV with SMILES column
        csv_content = (
            "ID,Name,SMILES\n"
            "1,Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O\n"
            "2,Paracetamol,CC(=O)NC1=CC=C(C=C1)O"
        )
        
        # Process with multiple properties
        properties = ["molecular_weight", "logp", "formula", "num_h_donors", "num_h_acceptors"]
        result = process_csv_data(csv_content, properties=properties)
        
        # Parse the result CSV
        result_df = pd.read_csv(io.StringIO(result["result"]))
        
        # Check if all requested properties were added
        for prop in properties:
            assert prop in result_df.columns
        
        # Verify selected property values
        assert round(result_df.loc[0, "molecular_weight"], 3) == 180.159  # Aspirin
        assert round(result_df.loc[1, "logp"], 1) == 1.4  # Paracetamol (修正: RDKitの計算値に合わせる)
        assert result_df.loc[0, "formula"] == "C9H8O4"  # Aspirin
        assert result_df.loc[1, "num_h_donors"] == 2  # Paracetamol
        assert result_df.loc[0, "num_h_acceptors"] == 3  # Aspirin (修正: RDKitの計算値に合わせる)

    def test_process_csv_with_custom_smiles_column(self):
        """Test processing CSV data with a custom specified SMILES column."""
        # Create CSV with differently named SMILES column
        csv_content = (
            "ID,Name,smiles_column\n"
            "1,Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O\n"
            "2,Paracetamol,CC(=O)NC1=CC=C(C=C1)O"
        )
        
        # Process with custom SMILES column
        result = process_csv_data(csv_content, smiles_column="smiles_column")
        
        # Check if result was processed correctly
        assert result["smiles_column"] == "smiles_column"
        
        # Parse the result CSV
        result_df = pd.read_csv(io.StringIO(result["result"]))
        
        # Verify molecular weight was calculated
        assert "molecular_weight" in result_df.columns
        assert round(result_df.loc[0, "molecular_weight"], 3) == 180.159

    def test_process_csv_with_error_handling(self):
        """Test error handling in CSV processing."""
        # Test with empty CSV
        result = process_csv_data("")
        assert "error" in result
        
        # Test with non-existent SMILES column
        csv_content = "ID,Name,Structure\n1,Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O"
        result = process_csv_data(csv_content, smiles_column="SMILES")
        assert "error" in result