"""
Module providing functionality for calculating molecular properties from molecular structures
"""
import logging
from typing import Dict, Union

# Logger configuration
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# RDKit import
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    rdkit_available = True
except ImportError:
    logger.warning("Unable to import RDKit module. Please verify it is installed on your system.")
    rdkit_available = False


def calculate_molecular_weight(smiles: str) -> float:
    """
    Calculate molecular weight from SMILES string

    Args:
        smiles: Molecular structure in SMILES notation

    Returns:
        float: Molecular weight. Returns NaN if SMILES string is invalid
    """
    if not rdkit_available:
        logger.error("Cannot calculate molecular weight because RDKit is not installed.")
        return float('nan')
        
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Invalid SMILES string: {smiles}")
            return float('nan')
        return Descriptors.MolWt(mol)
    except Exception as e:
        logger.error(f"Error occurred during molecular weight calculation: {str(e)}")
        return float('nan')


def calculate_properties(smiles: str) -> Dict[str, Union[float, str, None]]:
    """
    Get basic molecular properties from SMILES string

    Args:
        smiles: Molecular structure in SMILES notation

    Returns:
        Dict: Dictionary containing basic molecular properties
    """
    props = {
        "molecular_weight": float('nan'),
        "logp": float('nan'),
        "num_h_donors": float('nan'),
        "num_h_acceptors": float('nan'),
        "formula": None
    }
    
    if not rdkit_available:
        logger.error("Cannot calculate molecular properties because RDKit is not installed.")
        return props
        
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Invalid SMILES string: {smiles}")
            return props
            
        props["molecular_weight"] = Descriptors.MolWt(mol)
        props["logp"] = Descriptors.MolLogP(mol)
        props["num_h_donors"] = Descriptors.NumHDonors(mol)
        props["num_h_acceptors"] = Descriptors.NumHAcceptors(mol)
        props["formula"] = Chem.rdMolDescriptors.CalcMolFormula(mol)
    except Exception as e:
        logger.error(f"Error occurred during property calculation: {str(e)}")
        
    return props


# Adding Lipinski's Rule of Five check functionality for future expansion
def check_lipinski_rule(smiles: str) -> Dict[str, bool]:
    """
    Check Lipinski's Rule of Five for a SMILES string

    Args:
        smiles: Molecular structure in SMILES notation

    Returns:
        Dict: Results for each Lipinski's Rule condition
    """
    rules = {
        "molecular_weight_ok": False,  # ≤ 500
        "logp_ok": False,              # ≤ 5
        "h_donors_ok": False,          # ≤ 5
        "h_acceptors_ok": False,       # ≤ 10
        "all_rules_passed": False
    }
    
    if not rdkit_available:
        logger.error("Cannot check Lipinski's Rule because RDKit is not installed.")
        return rules
        
    try:
        props = calculate_properties(smiles)
        
        rules["molecular_weight_ok"] = props["molecular_weight"] <= 500
        rules["logp_ok"] = props["logp"] <= 5
        rules["h_donors_ok"] = props["num_h_donors"] <= 5
        rules["h_acceptors_ok"] = props["num_h_acceptors"] <= 10
        
        # Check if all rules are passed
        rules["all_rules_passed"] = all([
            rules["molecular_weight_ok"],
            rules["logp_ok"],
            rules["h_donors_ok"],
            rules["h_acceptors_ok"]
        ])
        
    except Exception as e:
        logger.error(f"Error occurred during Lipinski's Rule check: {str(e)}")
        
    return rules