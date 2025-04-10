"""
Module providing functionality for calculating molecular properties from molecular structures
"""
import logging
from typing import Dict, Union, List, Any, Optional, Tuple

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
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
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


def get_property_descriptions() -> Dict[str, Dict[str, str]]:
    """
    Get descriptions of all available molecular properties.
    
    Returns:
        Dict: Dictionary containing property names as keys and their descriptions as values
    """
    return {
        # Basic properties
        "molecular_weight": {
            "ja": "分子量",
            "en": "Molecular Weight",
            "description": "分子の平均分子量（各元素の平均原子量に基づく）。標準的な分子量を表す。",
            "module": "Descriptors.MolWt"
        },
        "exact_mol_wt": {
            "ja": "厳密分子量",
            "en": "Exact Molecular Weight",
            "description": "同位体組成を考慮した厳密な分子量（モノアイソトピック質量）。小数点まで厳密に算出された質量。",
            "module": "Descriptors.ExactMolWt"
        },
        "heavy_atom_mol_wt": {
            "ja": "重原子分子量",
            "en": "Heavy Atom Molecular Weight",
            "description": "水素を無視した分子の平均分子量。炭素など重原子のみの質量合計。",
            "module": "Descriptors.HeavyAtomMolWt"
        },
        "formula": {
            "ja": "分子式",
            "en": "Molecular Formula",
            "description": "分子の化学式。分子を構成する原子の種類と数を表す。",
            "module": "rdMolDescriptors.CalcMolFormula"
        },
        
        # Lipophilicity/Hydrophilicity
        "logp": {
            "ja": "オクタノール/水分配係数",
            "en": "LogP",
            "description": "分配係数LogP（1-オクタノールと水の間の分配係数の対数）。Wildman–Crippen法による推算値。疎水性の指標。",
            "module": "Descriptors.MolLogP"
        },
        "mol_mr": {
            "ja": "モル屈折率",
            "en": "Molar Refractivity",
            "description": "モル屈折率。Wildman–Crippen法による推算値で、分子の分極率に関連する指標。",
            "module": "Descriptors.MolMR"
        },
        "tpsa": {
            "ja": "トポロジカル極性表面積",
            "en": "Topological Polar Surface Area",
            "description": "分子中の極性原子（主にOとN、およびそれに結合したH）の表面積の合計。極性表面積が大きいほど細胞膜透過性は低下する傾向がある。",
            "module": "Descriptors.TPSA"
        },

        # Surface areas
        "labute_asa": {
            "ja": "Labute推定表面積",
            "en": "Labute's Approx. Surface Area",
            "description": "Labuteによる近似分子表面積。MOEソフトウェア由来のアルゴリズムで計算される分子表面積。",
            "module": "Descriptors.LabuteASA"
        },
        
        # H-bonds and atom counts
        "num_h_donors": {
            "ja": "水素結合ドナー数",
            "en": "#H-Bond Donors",
            "description": "分子内の水素結合供与体の個数。一般に–OHや–NH基の数（Lipinskiの定義に基づく）。",
            "module": "Descriptors.NumHDonors"
        },
        "num_h_acceptors": {
            "ja": "水素結合アクセプター数",
            "en": "#H-Bond Acceptors",
            "description": "分子内の水素結合受容体の個数。一般にOやN原子の数（Lipinskiの定義に基づく）。",
            "module": "Descriptors.NumHAcceptors"
        },
        "num_rotatable_bonds": {
            "ja": "回転可能結合数",
            "en": "#Rotatable Bonds",
            "description": "分子内の回転可能な単結合の数（末端の単結合や二重結合環境を除く）。分子の柔軟性の指標。",
            "module": "Descriptors.NumRotatableBonds"
        },
        "heavy_atom_count": {
            "ja": "重原子数",
            "en": "#Heavy Atoms",
            "description": "分子中の重原子（非水素原子）の数。分子サイズの指標の一つ。",
            "module": "Descriptors.HeavyAtomCount"
        },
        "num_hetero_atoms": {
            "ja": "ヘテロ原子数",
            "en": "#Heteroatoms",
            "description": "分子中のヘテロ原子（炭素以外の元素）の数。例えばN, O, Sなどの個数。",
            "module": "Descriptors.NumHeteroatoms"
        },
        "no_count": {
            "ja": "N/O原子数",
            "en": "#Nitrogen/Oxygen Atoms",
            "description": "分子中の窒素原子および酸素原子の合計数。LipinskiのH受容体数と近似的に対応。",
            "module": "Lipinski.NOCount"
        },
        "nhoh_count": {
            "ja": "OH/NH基数",
            "en": "#OH/NH Groups",
            "description": "分子中のヒドロキシ基(-OH)および一次・二次アミン基(-NH)の数。LipinskiのHドナー数と近い概念。",
            "module": "Lipinski.NHOHCount"
        },
        "num_valence_electrons": {
            "ja": "原子価電子数",
            "en": "#Valence Electrons",
            "description": "分子内の全原子の価電子の総数。分子全体の価電子の数で、電荷や元素組成の指標。",
            "module": "Descriptors.NumValenceElectrons"
        },
        
        # Ring information
        "num_aromatic_rings": {
            "ja": "芳香環数",
            "en": "#Aromatic Rings",
            "description": "分子内の芳香族環の数。ベンゼン環など芳香環構造の個数。",
            "module": "Descriptors.NumAromaticRings"
        },
        "num_aliphatic_rings": {
            "ja": "脂肪族環数",
            "en": "#Aliphatic Rings",
            "description": "分子内の脂肪族環の数。非芳香族の環構造（シクロアルカンなど）の個数。",
            "module": "Descriptors.NumAliphaticRings"
        },
        "num_saturated_rings": {
            "ja": "飽和環数",
            "en": "#Saturated Rings",
            "description": "分子内の飽和環の数。二重結合を含まない環（飽和脂肪族環）の個数。",
            "module": "Descriptors.NumSaturatedRings"
        },
        "num_aromatic_carbocycles": {
            "ja": "芳香族炭素環数",
            "en": "#Aromatic Carbocycles",
            "description": "芳香族炭素環（全原子が炭素の芳香環）の数。",
            "module": "Descriptors.NumAromaticCarbocycles"
        },
        "num_aromatic_heterocycles": {
            "ja": "芳香族複素環数",
            "en": "#Aromatic Heterocycles",
            "description": "芳香族複素環（炭素以外の原子を含む芳香環）の数。",
            "module": "Descriptors.NumAromaticHeterocycles"
        },
        "num_aliphatic_carbocycles": {
            "ja": "脂肪族炭素環数",
            "en": "#Aliphatic Carbocycles",
            "description": "脂肪族炭素環（炭素のみからなる非芳香環）の数。",
            "module": "Descriptors.NumAliphaticCarbocycles"
        },
        "num_aliphatic_heterocycles": {
            "ja": "脂肪族複素環数",
            "en": "#Aliphatic Heterocycles",
            "description": "脂肪族複素環（ヘテロ原子を含む非芳香環）の数。",
            "module": "Descriptors.NumAliphaticHeterocycles"
        },
        "num_saturated_carbocycles": {
            "ja": "飽和炭素環数",
            "en": "#Saturated Carbocycles",
            "description": "飽和炭素環（完全に単結合のみからなる炭素環）の数。",
            "module": "Descriptors.NumSaturatedCarbocycles"
        },
        "num_saturated_heterocycles": {
            "ja": "飽和複素環数",
            "en": "#Saturated Heterocycles",
            "description": "飽和複素環（完全に単結合のみからなる複素環）の数。",
            "module": "Descriptors.NumSaturatedHeterocycles"
        },
        "ring_count": {
            "ja": "環数（総数）",
            "en": "Ring Count",
            "description": "分子内の環構造の総数。全ての大小の環をカウントしたもの。",
            "module": "Descriptors.RingCount"
        },
        
        # Bond/functional group counts
        "num_amide_bonds": {
            "ja": "アミド結合数",
            "en": "#Amide Bonds",
            "description": "分子中のアミド結合（–C(=O)N–結合）の数。",
            "module": "Descriptors.NumAmideBonds"
        },
        "fraction_csp3": {
            "ja": "炭素sp³割合",
            "en": "Fraction of C sp³",
            "description": "炭素原子のうちsp³混成状態にあるものの割合。値域0〜1で、分子の飽和度を示す指標。",
            "module": "Descriptors.FractionCSP3"
        },
        
        # Molecular complexity
        "num_spiro_atoms": {
            "ja": "スピロ原子数",
            "en": "#Spiro Atoms",
            "description": "スピロ原子の数。2つの環が1つの原子を共有して融合した構造（スピロ環）における共有原子の個数。",
            "module": "Descriptors.NumSpiroAtoms"
        },
        "num_bridgehead_atoms": {
            "ja": "ブリッジヘッド原子数",
            "en": "#Bridgehead Atoms",
            "description": "ブリッジヘッド原子の数。複数の環が少なくとも2つの結合を共有して融合した構造（ブリッジヘッド）における共有原子の個数。",
            "module": "Descriptors.NumBridgeheadAtoms"
        },
        "num_stereo_centers": {
            "ja": "立体中心数",
            "en": "#Stereo Centers",
            "description": "分子内の立体中心（キラル中心）となっている原子の数。明確に(R/S等)指定された立体中心の個数。",
            "module": "Descriptors.NumAtomStereoCenters"
        },
        "num_unspecified_stereo_centers": {
            "ja": "未指定立体中心数",
            "en": "#Unspecified Stereo Centers",
            "description": "分子内の立体中心のうち立体配置が未定義のものの数（立体化学が指定されていないキラル原子の個数）。",
            "module": "Descriptors.NumUnspecifiedAtomStereoCenters"
        },
        
        # Graph indices
        "balaban_j": {
            "ja": "BalabanのJ指数",
            "en": "Balaban's J Index",
            "description": "Balabanが提唱した分子接続指数J。分子の接続グラフから計算されるトップロジカル指数で、分子の密結合度合いを表す。",
            "module": "GraphDescriptors.BalabanJ"
        },
        "bertz_ct": {
            "ja": "Bertzの複雑度指数",
            "en": "Bertz's Complexity Index",
            "description": "Bertzによる分子構造の複雑さ（Complexity）を定量化する指数。原子と結合の組み合わせの情報量に基づき分子の複雑度を評価する。",
            "module": "GraphDescriptors.BertzCT"
        },
        "ipc": {
            "ja": "Ipc情報指数",
            "en": "Ipc (Information Content)",
            "description": "分子グラフの情報量記述子Ipc。隣接行列の特性多項式係数から算出される情報量に基づく指標。",
            "module": "GraphDescriptors.Ipc"
        },
        "hall_kier_alpha": {
            "ja": "Hall–Kierのαパラメータ",
            "en": "Hall-Kier Alpha",
            "description": "Kierらによる分子の補正パラメータα。一般に分子のサイズや分枝に関する補正項で、κ形状指数の計算にも用いられる。",
            "module": "GraphDescriptors.HallKierAlpha"
        },

        # Kappa shape indices
        "kappa1": {
            "ja": "κ形状指数1",
            "en": "Kappa Shape Index 1",
            "description": "Kierの形状指数。与えられた原子数の中で最も線状な場合および最も分岐した場合を基準に、分子の形状を表す指数。Kappa1は分子の分岐度を表す。",
            "module": "GraphDescriptors.Kappa1"
        },
        "kappa2": {
            "ja": "κ形状指数2",
            "en": "Kappa Shape Index 2",
            "description": "Kierの形状指数。Kappa2は分子の空間的な広がり（面状）を表す。",
            "module": "GraphDescriptors.Kappa2"
        },
        "kappa3": {
            "ja": "κ形状指数3",
            "en": "Kappa Shape Index 3",
            "description": "Kierの形状指数。Kappa3は分子の空間的な広がり（立体的）を表す。",
            "module": "GraphDescriptors.Kappa3"
        },
        
        # Chi connectivity indices
        "chi0": {
            "ja": "χ0接続性指数",
            "en": "Chi0 Connectivity Index",
            "description": "KierとHallが提案した分子接続性指数。Chi0は0次の接続性指数。",
            "module": "GraphDescriptors.Chi0"
        },
        "chi1": {
            "ja": "χ1接続性指数",
            "en": "Chi1 Connectivity Index",
            "description": "KierとHallが提案した分子接続性指数。Chi1は1次の接続性指数。",
            "module": "GraphDescriptors.Chi1"
        },
        "chi0v": {
            "ja": "χ0v接続性指数（原子価考慮）",
            "en": "Chi0v Connectivity Index",
            "description": "KierとHallが提案した分子接続性指数。Chi0vは原子価を考慮した0次の接続性指数。",
            "module": "GraphDescriptors.Chi0v"
        },
        "chi1v": {
            "ja": "χ1v接続性指数（原子価考慮）",
            "en": "Chi1v Connectivity Index",
            "description": "KierとHallが提案した分子接続性指数。Chi1vは原子価を考慮した1次の接続性指数。",
            "module": "GraphDescriptors.Chi1v"
        },
        
        # QED drug-likeness
        "qed": {
            "ja": "QED薬剤様性スコア",
            "en": "QED Drug-Likeness Score",
            "description": "QED (Quantitative Estimation of Drug-likeness)。分子量、LogP、TPSA、Hドナー/アクセプター数、芳香環数、回転結合数など複数の性質を総合して0〜1の範囲で算出される薬剤様性指標。値が高いほどドラッグライクとされる。",
            "module": "QED.qed"
        },
    }


def get_properties_table(format: str = "markdown", language: str = "ja") -> str:
    """
    Generate a formatted table of all available molecular properties
    
    Args:
        format: Output format ('markdown', 'csv', 'json', 'text')
        language: Language for property names ('ja' for Japanese, 'en' for English)
    
    Returns:
        str: Formatted table of properties
    """
    descriptions = get_property_descriptions()
    
    if format == "json":
        import json
        # JSON形式で出力
        result = {}
        for prop_name, prop_info in descriptions.items():
            result[prop_name] = {
                "name": prop_info[language],
                "description": prop_info["description"],
                "module": prop_info["module"]
            }
        return json.dumps(result, indent=2, ensure_ascii=False)
    
    elif format == "csv":
        import csv
        import io
        # CSV形式で出力
        output = io.StringIO()
        fieldnames = ["property_id", "name", "description", "module"]
        writer = csv.DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()
        
        for prop_name, prop_info in descriptions.items():
            writer.writerow({
                "property_id": prop_name,
                "name": prop_info[language],
                "description": prop_info["description"],
                "module": prop_info["module"]
            })
        
        return output.getvalue()
    
    elif format == "text":
        # テキスト形式（単純なリスト）で出力
        lines = []
        for prop_name, prop_info in descriptions.items():
            lines.append(f"{prop_name} ({prop_info[language]}): {prop_info['description']}")
        return "\n".join(lines)
    
    else:  # デフォルト: markdown
        # Markdown形式でテーブルを出力
        lines = []
        lines.append("| プロパティID | 名称 | 説明 | モジュール |" if language == "ja" else 
                    "| Property ID | Name | Description | Module |")
        lines.append("|------------|------|------|----------|")
        
        for prop_name, prop_info in descriptions.items():
            lines.append(f"| {prop_name} | {prop_info[language]} | {prop_info['description']} | {prop_info['module']} |")
        
        return "\n".join(lines)


def get_property_categories() -> Dict[str, List[str]]:
    """
    Get molecular properties grouped by categories
    
    Returns:
        Dict: Dictionary with categories as keys and lists of property IDs as values
    """
    return {
        "basic": ["molecular_weight", "exact_mol_wt", "heavy_atom_mol_wt", "formula"],
        "lipophilicity": ["logp", "mol_mr", "tpsa"],
        "surface_area": ["labute_asa"],
        "h_bonds_and_atoms": ["num_h_donors", "num_h_acceptors", "num_rotatable_bonds", "heavy_atom_count", 
                            "num_hetero_atoms", "no_count", "nhoh_count", "num_valence_electrons"],
        "rings": ["num_aromatic_rings", "num_aliphatic_rings", "num_saturated_rings",
                "num_aromatic_carbocycles", "num_aromatic_heterocycles", 
                "num_aliphatic_carbocycles", "num_aliphatic_heterocycles",
                "num_saturated_carbocycles", "num_saturated_heterocycles", "ring_count"],
        "bonds_and_groups": ["num_amide_bonds", "fraction_csp3"],
        "complexity": ["num_spiro_atoms", "num_bridgehead_atoms", "num_stereo_centers", 
                     "num_unspecified_stereo_centers"],
        "graph_indices": ["balaban_j", "bertz_ct", "ipc", "hall_kier_alpha", 
                        "kappa1", "kappa2", "kappa3", "chi0", "chi1", "chi0v", "chi1v"],
        "drug_likeness": ["qed"]
    }


def get_available_properties() -> List[str]:
    """
    Get a list of all available property names that can be calculated
    
    Returns:
        List[str]: List of all property names that can be calculated
    """
    # サンプルのSMILES文字列を使用して全プロパティ名を取得
    descriptions = get_property_descriptions()
    return list(descriptions.keys())


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


def check_pains_filter(smiles: str) -> Dict[str, Any]:
    """
    Check if a molecule contains PAINS (Pan Assay Interference Compounds) substructures
    
    Args:
        smiles: Molecular structure in SMILES notation
        
    Returns:
        Dict: Results of PAINS filter check with details about matching patterns
    """
    result = {
        "pains_free": True,        # PAINSアラートパターンがない
        "num_alerts": 0,           # 検出されたアラート数
        "matching_alerts": [],     # マッチしたアラートの詳細
    }
    
    if not rdkit_available:
        logger.error("Cannot check PAINS filter because RDKit is not installed.")
        return result
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Invalid SMILES string: {smiles}")
            return result
        
        # PAINS用のフィルターカタログを作成
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
        catalog = FilterCatalog(params)
        
        # 分子にPAINSパターンがあるかチェック
        if catalog.HasMatch(mol):
            # マッチしたエントリを取得
            matches = catalog.GetMatches(mol)
            result["pains_free"] = False
            result["num_alerts"] = len(matches)
            
            # マッチしたアラートの詳細情報を取得
            for match in matches:
                match_dict = {
                    "description": match.GetDescription(),
                    "smarts": match.GetSmarts() if hasattr(match, "GetSmarts") else None
                }
                result["matching_alerts"].append(match_dict)
                
    except Exception as e:
        logger.error(f"Error occurred during PAINS filter check: {str(e)}")
    
    return result


def check_egan_filter(smiles: str) -> Dict[str, bool]:
    """
    Check Egan filter for drug-likeness (good oral absorption)
    
    Args:
        smiles: Molecular structure in SMILES notation
        
    Returns:
        Dict: Results for each Egan filter condition
    """
    rules = {
        "logp_ok": False,         # ≤ 5.88
        "tpsa_ok": False,         # ≤ 131.6
        "all_rules_passed": False
    }
    
    if not rdkit_available:
        logger.error("Cannot check Egan filter because RDKit is not installed.")
        return rules
    
    try:
        # 必要なプロパティを計算
        props = calculate_selected_properties(smiles, ["logp", "tpsa"])
        
        logp = props.get("logp", float('inf'))
        tpsa = props.get("tpsa", float('inf'))
        
        rules["logp_ok"] = logp <= 5.88
        rules["tpsa_ok"] = tpsa <= 131.6
        
        # 全てのルールをパスしたかチェック
        rules["all_rules_passed"] = all([
            rules["logp_ok"],
            rules["tpsa_ok"]
        ])
        
    except Exception as e:
        logger.error(f"Error occurred during Egan filter check: {str(e)}")
    
    return rules


def check_muegge_filter(smiles: str) -> Dict[str, bool]:
    """
    Check Muegge filter for drug-likeness
    
    Args:
        smiles: Molecular structure in SMILES notation
        
    Returns:
        Dict: Results for each Muegge filter condition
    """
    rules = {
        "molecular_weight_ok": False,     # 200-600
        "logp_ok": False,                 # -2-5
        "tpsa_ok": False,                 # ≤ 150
        "ring_count_ok": False,           # ≤ 7
        "h_acceptors_ok": False,          # ≤ 10
        "h_donors_ok": False,             # ≤ 5
        "rotatable_bonds_ok": False,      # < 15
        "all_rules_passed": False
    }
    
    if not rdkit_available:
        logger.error("Cannot check Muegge filter because RDKit is not installed.")
        return rules
    
    try:
        # 必要なプロパティを計算
        props = calculate_selected_properties(smiles, [
            "molecular_weight", "logp", "tpsa", "ring_count",
            "num_h_acceptors", "num_h_donors", "num_rotatable_bonds"
        ])
        
        mw = props.get("molecular_weight", 0)
        logp = props.get("logp", 0)
        tpsa = props.get("tpsa", float('inf'))
        ring_count = props.get("ring_count", float('inf'))
        h_acceptors = props.get("num_h_acceptors", float('inf'))
        h_donors = props.get("num_h_donors", float('inf'))
        rotatable_bonds = props.get("num_rotatable_bonds", float('inf'))
        
        rules["molecular_weight_ok"] = 200 <= mw <= 600
        rules["logp_ok"] = -2 <= logp <= 5
        rules["tpsa_ok"] = tpsa <= 150
        rules["ring_count_ok"] = ring_count <= 7
        rules["h_acceptors_ok"] = h_acceptors <= 10
        rules["h_donors_ok"] = h_donors <= 5
        rules["rotatable_bonds_ok"] = rotatable_bonds < 15
        
        # 全てのルールをパスしたかチェック
        rules["all_rules_passed"] = all([
            rules["molecular_weight_ok"],
            rules["logp_ok"],
            rules["tpsa_ok"],
            rules["ring_count_ok"],
            rules["h_acceptors_ok"],
            rules["h_donors_ok"],
            rules["rotatable_bonds_ok"]
        ])
        
    except Exception as e:
        logger.error(f"Error occurred during Muegge filter check: {str(e)}")
    
    return rules


def check_all_druglikeness_filters(smiles: str) -> Dict[str, Any]:
    """
    Run all available druglikeness filters on a molecule
    
    Args:
        smiles: Molecular structure in SMILES notation
        
    Returns:
        Dict: Results of all drug-likeness filter checks
    """
    result = {
        "lipinski": None,
        "veber": None,
        "ghose": None,
        "egan": None,
        "muegge": None,
        "pains": None,
        "all_filters_passed": False
    }
    
    if not rdkit_available:
        logger.error("Cannot check drug-likeness filters because RDKit is not installed.")
        return result
    
    try:
        # すべてのフィルターを実行
        result["lipinski"] = check_lipinski_rule(smiles)
        result["veber"] = check_veber_rules(smiles)
        result["ghose"] = check_ghose_filter(smiles)
        result["egan"] = check_egan_filter(smiles)
        result["muegge"] = check_muegge_filter(smiles)
        result["pains"] = check_pains_filter(smiles)
        
        # 全てのフィルターをパスしたかチェック
        result["all_filters_passed"] = (
            result["lipinski"]["all_rules_passed"] and
            result["veber"]["all_rules_passed"] and
            result["ghose"]["all_rules_passed"] and
            result["egan"]["all_rules_passed"] and
            result["muegge"]["all_rules_passed"] and
            result["pains"]["pains_free"]
        )
        
    except Exception as e:
        logger.error(f"Error occurred during drug-likeness filter checks: {str(e)}")
    
    return result


def get_druglikeness_filters_summary(smiles: str) -> Dict[str, bool]:
    """
    Get a simple summary of all drug-likeness filters for a molecule
    
    Args:
        smiles: Molecular structure in SMILES notation
        
    Returns:
        Dict: Simple summary of all drug-likeness filters with just pass/fail results
    """
    result = {
        "lipinski_pass": False,
        "veber_pass": False,
        "ghose_pass": False,
        "egan_pass": False, 
        "muegge_pass": False,
        "pains_free": False,
        "all_filters_passed": False
    }
    
    if not rdkit_available:
        logger.error("Cannot check drug-likeness filters because RDKit is not installed.")
        return result
    
    try:
        all_results = check_all_druglikeness_filters(smiles)
        
        result["lipinski_pass"] = all_results["lipinski"]["all_rules_passed"]
        result["veber_pass"] = all_results["veber"]["all_rules_passed"]
        result["ghose_pass"] = all_results["ghose"]["all_rules_passed"]
        result["egan_pass"] = all_results["egan"]["all_rules_passed"]
        result["muegge_pass"] = all_results["muegge"]["all_rules_passed"]
        result["pains_free"] = all_results["pains"]["pains_free"]
        result["all_filters_passed"] = all_results["all_filters_passed"]
        
    except Exception as e:
        logger.error(f"Error occurred during drug-likeness summary generation: {str(e)}")
    
    return result