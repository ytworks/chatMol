"""
分子構造から物性値を計算する機能を提供するモジュール
"""
import logging
from typing import Dict, Union

# ロガーの設定
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# RDKitのインポート
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    rdkit_available = True
except ImportError:
    logger.warning("RDKitモジュールをインポートできません。システムにインストールされているか確認してください。")
    rdkit_available = False


def calculate_molecular_weight(smiles: str) -> float:
    """
    SMILES文字列から分子量を計算する

    Args:
        smiles: SMILES表記の分子構造

    Returns:
        float: 分子量。SMILES文字列が無効の場合はNaN
    """
    if not rdkit_available:
        logger.error("RDKitがインストールされていないため、分子量を計算できません。")
        return float('nan')
        
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"無効なSMILES文字列: {smiles}")
            return float('nan')
        return Descriptors.MolWt(mol)
    except Exception as e:
        logger.error(f"分子量計算中にエラーが発生しました: {str(e)}")
        return float('nan')


def calculate_properties(smiles: str) -> Dict[str, Union[float, str, None]]:
    """
    SMILES文字列から分子の基本的な特性を取得

    Args:
        smiles: SMILES表記の分子構造

    Returns:
        Dict: 分子の基本特性を含む辞書
    """
    props = {
        "molecular_weight": float('nan'),
        "logp": float('nan'),
        "num_h_donors": float('nan'),
        "num_h_acceptors": float('nan'),
        "formula": None
    }
    
    if not rdkit_available:
        logger.error("RDKitがインストールされていないため、物性値を計算できません。")
        return props
        
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"無効なSMILES文字列: {smiles}")
            return props
            
        props["molecular_weight"] = Descriptors.MolWt(mol)
        props["logp"] = Descriptors.MolLogP(mol)
        props["num_h_donors"] = Descriptors.NumHDonors(mol)
        props["num_h_acceptors"] = Descriptors.NumHAcceptors(mol)
        props["formula"] = Chem.rdMolDescriptors.CalcMolFormula(mol)
    except Exception as e:
        logger.error(f"物性値計算中にエラーが発生しました: {str(e)}")
        
    return props


# 将来的な拡張のために、Lipinski's Rule of Fiveチェック機能を追加
def check_lipinski_rule(smiles: str) -> Dict[str, bool]:
    """
    SMILES文字列に対してLipinski's Rule of Fiveをチェックする

    Args:
        smiles: SMILES表記の分子構造

    Returns:
        Dict: Lipinski's Ruleの各条件に対する結果
    """
    rules = {
        "molecular_weight_ok": False,  # 500以下
        "logp_ok": False,              # 5以下
        "h_donors_ok": False,          # 5以下
        "h_acceptors_ok": False,       # 10以下
        "all_rules_passed": False
    }
    
    if not rdkit_available:
        logger.error("RDKitがインストールされていないため、Lipinski's Ruleをチェックできません。")
        return rules
        
    try:
        props = calculate_properties(smiles)
        
        rules["molecular_weight_ok"] = props["molecular_weight"] <= 500
        rules["logp_ok"] = props["logp"] <= 5
        rules["h_donors_ok"] = props["num_h_donors"] <= 5
        rules["h_acceptors_ok"] = props["num_h_acceptors"] <= 10
        
        # すべてのルールをパスしているかチェック
        rules["all_rules_passed"] = all([
            rules["molecular_weight_ok"],
            rules["logp_ok"],
            rules["h_donors_ok"],
            rules["h_acceptors_ok"]
        ])
        
    except Exception as e:
        logger.error(f"Lipinski's Ruleチェック中にエラーが発生しました: {str(e)}")
        
    return rules