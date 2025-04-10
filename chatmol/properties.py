"""
Module providing functionality for calculating molecular properties from molecular structures
"""
import logging
from typing import Dict, Union, List, Any, Optional

# Logger configuration
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# RDKit import
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, EState, QED, AllChem
    from rdkit.Chem import Crippen, MolSurf, rdMolDescriptors, Fragments
    from rdkit.Chem.EState import EState_VSA
    from rdkit.Chem import GraphDescriptors
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
        props["formula"] = rdMolDescriptors.CalcMolFormula(mol)
    except Exception as e:
        logger.error(f"Error occurred during property calculation: {str(e)}")
        
    return props


def calculate_all_properties(smiles: str) -> Dict[str, Union[float, str, None]]:
    """
    Calculate comprehensive molecular properties from SMILES string
    
    Args:
        smiles: Molecular structure in SMILES notation
        
    Returns:
        Dict: Dictionary containing comprehensive molecular properties
    """
    # Initialize with basic properties and NaN values
    props = {
        # Basic properties
        "molecular_weight": float('nan'),           # 分子量
        "exact_mol_wt": float('nan'),              # 厳密分子量
        "heavy_atom_mol_wt": float('nan'),         # 重原子分子量
        "formula": None,                           # 分子式
        
        # Lipophilicity/Hydrophilicity
        "logp": float('nan'),                      # LogP
        "mol_mr": float('nan'),                    # モル屈折率
        "tpsa": float('nan'),                      # トポロジカル極性表面積

        # Surface areas
        "labute_asa": float('nan'),                # Labute推定表面積
        
        # H-bonds and atom counts
        "num_h_donors": float('nan'),              # 水素結合ドナー数
        "num_h_acceptors": float('nan'),           # 水素結合アクセプター数
        "num_rotatable_bonds": float('nan'),       # 回転可能結合数
        "heavy_atom_count": float('nan'),          # 重原子数
        "num_hetero_atoms": float('nan'),          # ヘテロ原子数
        "no_count": float('nan'),                  # N/O原子数
        "nhoh_count": float('nan'),                # OH/NH基数
        "num_valence_electrons": float('nan'),     # 原子価電子数
        
        # Ring information
        "num_aromatic_rings": float('nan'),        # 芳香環数
        "num_aliphatic_rings": float('nan'),       # 脂肪族環数
        "num_saturated_rings": float('nan'),       # 飽和環数
        "num_aromatic_carbocycles": float('nan'),  # 芳香族炭素環数
        "num_aromatic_heterocycles": float('nan'), # 芳香族複素環数
        "num_aliphatic_carbocycles": float('nan'), # 脂肪族炭素環数
        "num_aliphatic_heterocycles": float('nan'),# 脂肪族複素環数
        "num_saturated_carbocycles": float('nan'), # 飽和炭素環数
        "num_saturated_heterocycles": float('nan'),# 飽和複素環数
        "ring_count": float('nan'),                # 環数（総数）
        
        # Bond/functional group counts
        "num_amide_bonds": float('nan'),           # アミド結合数
        "fraction_csp3": float('nan'),             # 炭素sp³割合
        
        # Molecular complexity
        "num_spiro_atoms": float('nan'),           # スピロ原子数
        "num_bridgehead_atoms": float('nan'),      # ブリッジヘッド原子数
        "num_stereo_centers": float('nan'),        # 立体中心数
        "num_unspecified_stereo_centers": float('nan'), # 未指定立体中心数
        
        # Charge related
        "max_partial_charge": float('nan'),        # 最大部分電荷
        "min_partial_charge": float('nan'),        # 最小部分電荷
        "max_abs_partial_charge": float('nan'),    # 最大絶対値部分電荷
        "min_abs_partial_charge": float('nan'),    # 最小絶対値部分電荷
        
        # EState indices
        "max_estate_index": float('nan'),          # 最大EState指数
        "min_estate_index": float('nan'),          # 最小EState指数
        "max_abs_estate_index": float('nan'),      # 最大絶対値EState指数
        "min_abs_estate_index": float('nan'),      # 最小絶対値EState指数

        # Graph indices
        "balaban_j": float('nan'),                 # BalabanのJ指数
        "bertz_ct": float('nan'),                  # Bertzの複雑度指数
        "ipc": float('nan'),                       # Ipc情報指数
        "hall_kier_alpha": float('nan'),           # Hall–Kierのαパラメータ
        
        # Kappa shape indices
        "kappa1": float('nan'),                    # κ形状指数1
        "kappa2": float('nan'),                    # κ形状指数2 
        "kappa3": float('nan'),                    # κ形状指数3
        
        # Chi connectivity indices
        "chi0": float('nan'),                      # χ0接続性指数
        "chi1": float('nan'),                      # χ1接続性指数
        "chi0v": float('nan'),                     # χ0v接続性指数（原子価考慮）
        "chi1v": float('nan'),                     # χ1v接続性指数（原子価考慮）
        
        # QED drug-likeness
        "qed": float('nan'),                       # QED薬剤様性スコア
    }
    
    # フラグメントカウント系のプロパティを初期化
    # あえて辞書に含めない（実行時に動的に処理する）
    
    if not rdkit_available:
        logger.error("Cannot calculate molecular properties because RDKit is not installed.")
        return props
        
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Invalid SMILES string: {smiles}")
            return props
            
        # 基本的な分子特性
        props["molecular_weight"] = Descriptors.MolWt(mol)
        props["exact_mol_wt"] = Descriptors.ExactMolWt(mol)
        props["heavy_atom_mol_wt"] = Descriptors.HeavyAtomMolWt(mol)
        props["formula"] = rdMolDescriptors.CalcMolFormula(mol)
        
        # 親油性/親水性
        props["logp"] = Descriptors.MolLogP(mol)
        props["mol_mr"] = Descriptors.MolMR(mol)
        props["tpsa"] = Descriptors.TPSA(mol)
        
        # 表面積
        props["labute_asa"] = Descriptors.LabuteASA(mol)
        
        # 水素結合と原子数
        props["num_h_donors"] = Descriptors.NumHDonors(mol)
        props["num_h_acceptors"] = Descriptors.NumHAcceptors(mol)
        props["num_rotatable_bonds"] = Descriptors.NumRotatableBonds(mol)
        props["heavy_atom_count"] = Descriptors.HeavyAtomCount(mol)
        props["num_hetero_atoms"] = Descriptors.NumHeteroatoms(mol)
        props["no_count"] = Lipinski.NOCount(mol)
        props["nhoh_count"] = Lipinski.NHOHCount(mol)
        props["num_valence_electrons"] = Descriptors.NumValenceElectrons(mol)
        
        # 環構造情報
        props["num_aromatic_rings"] = Descriptors.NumAromaticRings(mol)
        props["num_aliphatic_rings"] = Descriptors.NumAliphaticRings(mol)
        props["num_saturated_rings"] = Descriptors.NumSaturatedRings(mol)
        props["num_aromatic_carbocycles"] = Descriptors.NumAromaticCarbocycles(mol)
        props["num_aromatic_heterocycles"] = Descriptors.NumAromaticHeterocycles(mol)
        props["num_aliphatic_carbocycles"] = Descriptors.NumAliphaticCarbocycles(mol)
        props["num_aliphatic_heterocycles"] = Descriptors.NumAliphaticHeterocycles(mol)
        props["num_saturated_carbocycles"] = Descriptors.NumSaturatedCarbocycles(mol)
        props["num_saturated_heterocycles"] = Descriptors.NumSaturatedHeterocycles(mol)
        props["ring_count"] = Descriptors.RingCount(mol)
        
        # 結合/官能基カウント
        props["num_amide_bonds"] = Descriptors.NumAmideBonds(mol)
        props["fraction_csp3"] = Descriptors.FractionCSP3(mol)
        
        # 分子複雑性
        props["num_spiro_atoms"] = Descriptors.NumSpiroAtoms(mol)
        props["num_bridgehead_atoms"] = Descriptors.NumBridgeheadAtoms(mol)
        props["num_stereo_centers"] = Descriptors.NumAtomStereoCenters(mol)
        props["num_unspecified_stereo_centers"] = Descriptors.NumUnspecifiedAtomStereoCenters(mol)
        
        # 電荷関連（部分電荷計算要求がある場合）
        try:
            # 部分電荷計算の試み（失敗しても処理継続）
            AllChem.ComputeGasteigerCharges(mol)
            charges = [float(atom.GetProp('_GasteigerCharge')) for atom in mol.GetAtoms()]
            if charges:
                props["max_partial_charge"] = max(charges)
                props["min_partial_charge"] = min(charges)
                abs_charges = [abs(c) for c in charges]
                props["max_abs_partial_charge"] = max(abs_charges)
                props["min_abs_partial_charge"] = min(abs_charges)
        except Exception as e:
            logger.warning(f"Failed to compute partial charges: {str(e)}")

        # EState指数関連
        try:
            estate_indices = EState.EStateIndices(mol)
            if estate_indices:
                props["max_estate_index"] = max(estate_indices)
                props["min_estate_index"] = min(estate_indices)
                abs_estate = [abs(i) for i in estate_indices]
                props["max_abs_estate_index"] = max(abs_estate)
                props["min_abs_estate_index"] = min(abs_estate)
        except Exception as e:
            logger.warning(f"Failed to compute EState indices: {str(e)}")
            
        # グラフ指数
        props["balaban_j"] = GraphDescriptors.BalabanJ(mol)
        props["bertz_ct"] = GraphDescriptors.BertzCT(mol)
        props["ipc"] = GraphDescriptors.Ipc(mol)
        props["hall_kier_alpha"] = GraphDescriptors.HallKierAlpha(mol)

        # κ形状指数
        props["kappa1"] = GraphDescriptors.Kappa1(mol)
        props["kappa2"] = GraphDescriptors.Kappa2(mol)
        props["kappa3"] = GraphDescriptors.Kappa3(mol)
        
        # χ接続性指数
        props["chi0"] = GraphDescriptors.Chi0(mol)
        props["chi1"] = GraphDescriptors.Chi1(mol)
        props["chi0v"] = GraphDescriptors.Chi0v(mol)
        props["chi1v"] = GraphDescriptors.Chi1v(mol)

        # QED薬剤様性
        try:
            props["qed"] = QED.qed(mol)
        except Exception as e:
            logger.warning(f"Failed to compute QED: {str(e)}")
        
        # フラグメント解析（官能基カウント）
        fragment_props = {}
        for name in dir(Fragments):
            if name.startswith('fr_'):
                try:
                    func = getattr(Fragments, name)
                    if callable(func):
                        fragment_props[name] = func(mol)
                except Exception as e:
                    logger.debug(f"Failed to compute {name}: {str(e)}")
        
        # フラグメント解析結果をpropsに追加
        props.update(fragment_props)
        
    except Exception as e:
        logger.error(f"Error occurred during comprehensive property calculation: {str(e)}")
        
    return props


def calculate_selected_properties(smiles: str, properties_list: List[str]) -> Dict[str, Any]:
    """
    Calculate only selected molecular properties from SMILES string
    
    Args:
        smiles: Molecular structure in SMILES notation
        properties_list: List of property names to calculate
        
    Returns:
        Dict: Dictionary containing only the requested molecular properties
    """
    # まず全ての物性値を計算
    all_props = calculate_all_properties(smiles)
    
    # 要求されたプロパティのみを抽出
    selected_props = {prop: all_props.get(prop, None) for prop in properties_list if prop in all_props}
    
    # 要求されたが利用できないプロパティをログ出力
    missing_props = [prop for prop in properties_list if prop not in all_props]
    if missing_props:
        logger.warning(f"The following requested properties are not available: {', '.join(missing_props)}")
    
    return selected_props


def get_available_properties() -> List[str]:
    """
    Get a list of all available property names that can be calculated
    
    Returns:
        List[str]: List of all property names that can be calculated
    """
    # サンプルのSMILES文字列を使用して全プロパティ名を取得
    sample_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # アスピリン
    try:
        sample_props = calculate_all_properties(sample_smiles)
        return list(sample_props.keys())
    except Exception:
        # 基本的なプロパティリストを返す
        return [
            "molecular_weight", "logp", "num_h_donors", "num_h_acceptors", "formula",
            "tpsa", "heavy_atom_count", "num_rotatable_bonds"
        ]


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


def check_veber_rules(smiles: str) -> Dict[str, bool]:
    """
    Check Veber's Rules for oral bioavailability
    
    Args:
        smiles: Molecular structure in SMILES notation
        
    Returns:
        Dict: Results for each Veber's Rule condition
    """
    rules = {
        "rotatable_bonds_ok": False,  # ≤ 10
        "tpsa_ok": False,             # ≤ 140
        "all_rules_passed": False
    }
    
    if not rdkit_available:
        logger.error("Cannot check Veber's Rules because RDKit is not installed.")
        return rules
    
    try:
        # 必要なプロパティを計算
        props = calculate_selected_properties(smiles, ["num_rotatable_bonds", "tpsa"])
        
        rules["rotatable_bonds_ok"] = props.get("num_rotatable_bonds", float('inf')) <= 10
        rules["tpsa_ok"] = props.get("tpsa", float('inf')) <= 140
        
        # 全てのルールをパスしたかチェック
        rules["all_rules_passed"] = all([
            rules["rotatable_bonds_ok"],
            rules["tpsa_ok"]
        ])
        
    except Exception as e:
        logger.error(f"Error occurred during Veber's Rules check: {str(e)}")
    
    return rules


def check_ghose_filter(smiles: str) -> Dict[str, bool]:
    """
    Check Ghose filter for drug-likeness
    
    Args:
        smiles: Molecular structure in SMILES notation
        
    Returns:
        Dict: Results for each Ghose filter condition
    """
    rules = {
        "molecular_weight_ok": False,  # 160-480
        "logp_ok": False,              # -0.4-5.6
        "atom_count_ok": False,        # 20-70 atoms
        "molar_refractivity_ok": False, # 40-130
        "all_rules_passed": False
    }
    
    if not rdkit_available:
        logger.error("Cannot check Ghose filter because RDKit is not installed.")
        return rules
    
    try:
        # 必要なプロパティを計算
        props = calculate_selected_properties(smiles, [
            "molecular_weight", "logp", "heavy_atom_count", "mol_mr"
        ])
        
        mw = props.get("molecular_weight", 0)
        logp = props.get("logp", 0)
        atom_count = props.get("heavy_atom_count", 0)
        mol_mr = props.get("mol_mr", 0)
        
        rules["molecular_weight_ok"] = 160 <= mw <= 480
        rules["logp_ok"] = -0.4 <= logp <= 5.6
        rules["atom_count_ok"] = 20 <= atom_count <= 70
        rules["molar_refractivity_ok"] = 40 <= mol_mr <= 130
        
        # 全てのルールをパスしたかチェック
        rules["all_rules_passed"] = all([
            rules["molecular_weight_ok"],
            rules["logp_ok"],
            rules["atom_count_ok"],
            rules["molar_refractivity_ok"]
        ])
        
    except Exception as e:
        logger.error(f"Error occurred during Ghose filter check: {str(e)}")
    
    return rules