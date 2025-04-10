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
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
    rdkit_available = True
except ImportError:
    logger.warning("Unable to import RDKit module. Please verify it is installed on your system.")
    rdkit_available = False


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


def get_available_properties() -> List[str]:
    """
    Get a list of all available property names that can be calculated
    
    Returns:
        List[str]: List of all property names that can be calculated
    """
    # サンプルのSMILES文字列を使用して全プロパティ名を取得
    descriptions = get_property_descriptions()
    return list(descriptions.keys())


def get_feature_descriptions() -> Dict[str, Dict[str, str]]:
    """
    すべての分子特性（プロパティとフィルター）の説明を取得する
    
    Returns:
        Dict: 分子特性とその説明の辞書
    """
    # プロパティの説明を取得
    feature_descriptions = get_property_descriptions()
    
    # 元のMOLECULAR_FILTERSの情報を新しい形式で統合
    for filter_name, filter_info in MOLECULAR_FILTERS.items():
        # 既存のプロパティ説明と同じ形式に変換
        feature_descriptions[filter_name] = {
            "ja": filter_name,  # 日本語名がない場合はキー名をそのまま使用
            "en": filter_name,  # 英語名がない場合はキー名をそのまま使用
            "description": filter_info.get("description", ""),
            "is_filter": True,  # フィルターであることを示すフラグを追加
            "result_key": filter_info.get("result_key", "")  # フィルター固有の情報を保持
        }
    
    return feature_descriptions


# フィルターとその結果キーのマッピング
MOLECULAR_FILTERS = {
    "lipinski": {
        "result_key": "all_rules_passed",
        "description": "Lipinski's Rule of Five (MW≤500, LogP≤5, HBD≤5, HBA≤10)"
    },
    "veber": {
        "result_key": "all_rules_passed",
        "description": "Veber's rules (TPSA≤140, RotBonds≤10)"
    },
    "ghose": {
        "result_key": "all_rules_passed",
        "description": "Ghose filter (160≤MW≤480, -0.4≤LogP≤5.6, 20≤atoms≤70, 40≤MR≤130)"
    },
    "egan": {
        "result_key": "all_rules_passed",
        "description": "Egan filter (LogP≤5.88, TPSA≤131.6)"
    },
    "muegge": {
        "result_key": "all_rules_passed",
        "description": "Muegge filter (200≤MW≤600, -2≤LogP≤5, TPSA≤150, rings≤7, HBA≤10, HBD≤5, RotBonds<15)"
    },
    "pains": {
        "result_key": "pains_free",
        "description": "PAINS filter (screens for pan-assay interference compounds)"
    }
}


def calculate_molecular_features(
    smiles: str, 
    use_rdkit_mol: bool = False
) -> Dict[str, Any]:
    """
    分子の全ての物性値とフィルター判定を計算しフラットな辞書で返す統合関数
    
    Args:
        smiles: 分子構造をSMILES表記で
        use_rdkit_mol: Trueの場合、RDKit分子オブジェクトも返す（再利用のため）
    
    Returns:
        Dict: 計算された全ての物性値とフィルター結果をフラットに含む辞書
    """
    # 結果辞書を初期化（フラットな構造）
    result = {
        "smiles": smiles
    }
    
    # RDKitが利用可能かチェック
    if not rdkit_available:
        logger.error("RDKitがインストールされていないため、分子特性を計算できません")
        if use_rdkit_mol:
            result["mol"] = None
        return result
    
    try:
        # SMILES文字列からRDKit分子オブジェクトを作成
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"無効なSMILES文字列です: {smiles}")
            if use_rdkit_mol:
                result["mol"] = None
            return result
        
        # 分子オブジェクトを保存（オプション）
        if use_rdkit_mol:
            result["mol"] = mol
        
        # 基本的な分子特性
        result["molecular_weight"] = Descriptors.MolWt(mol)
        result["exact_mol_wt"] = Descriptors.ExactMolWt(mol)
        result["heavy_atom_mol_wt"] = Descriptors.HeavyAtomMolWt(mol)
        result["formula"] = rdMolDescriptors.CalcMolFormula(mol)
        
        # 親油性/親水性
        result["logp"] = Descriptors.MolLogP(mol)
        result["mol_mr"] = Descriptors.MolMR(mol)
        result["tpsa"] = Descriptors.TPSA(mol)
        
        # 表面積
        result["labute_asa"] = Descriptors.LabuteASA(mol)
        
        # 水素結合と原子数
        result["num_h_donors"] = Descriptors.NumHDonors(mol)
        result["num_h_acceptors"] = Descriptors.NumHAcceptors(mol)
        result["num_rotatable_bonds"] = Descriptors.NumRotatableBonds(mol)
        result["heavy_atom_count"] = Descriptors.HeavyAtomCount(mol)
        result["num_hetero_atoms"] = Descriptors.NumHeteroatoms(mol)
        result["no_count"] = Lipinski.NOCount(mol)
        result["nhoh_count"] = Lipinski.NHOHCount(mol)
        result["num_valence_electrons"] = Descriptors.NumValenceElectrons(mol)
        
        # 環構造情報
        result["num_aromatic_rings"] = Descriptors.NumAromaticRings(mol)
        result["num_aliphatic_rings"] = Descriptors.NumAliphaticRings(mol)
        result["num_saturated_rings"] = Descriptors.NumSaturatedRings(mol)
        result["num_aromatic_carbocycles"] = Descriptors.NumAromaticCarbocycles(mol)
        result["num_aromatic_heterocycles"] = Descriptors.NumAromaticHeterocycles(mol)
        result["num_aliphatic_carbocycles"] = Descriptors.NumAliphaticCarbocycles(mol)
        result["num_aliphatic_heterocycles"] = Descriptors.NumAliphaticHeterocycles(mol)
        result["num_saturated_carbocycles"] = Descriptors.NumSaturatedCarbocycles(mol)
        result["num_saturated_heterocycles"] = Descriptors.NumSaturatedHeterocycles(mol)
        result["ring_count"] = Descriptors.RingCount(mol)
        
        # 結合/官能基カウント
        result["num_amide_bonds"] = Descriptors.NumAmideBonds(mol)
        result["fraction_csp3"] = Descriptors.FractionCSP3(mol)
        
        # 分子複雑性
        result["num_spiro_atoms"] = Descriptors.NumSpiroAtoms(mol)
        result["num_bridgehead_atoms"] = Descriptors.NumBridgeheadAtoms(mol)
        result["num_stereo_centers"] = Descriptors.NumAtomStereoCenters(mol)
        result["num_unspecified_stereo_centers"] = Descriptors.NumUnspecifiedAtomStereoCenters(mol)
        
        # 電荷関連の物性値
        try:
            AllChem.ComputeGasteigerCharges(mol)
            charges = [float(atom.GetProp('_GasteigerCharge')) for atom in mol.GetAtoms()]
            if charges:
                result["max_partial_charge"] = max(charges)
                result["min_partial_charge"] = min(charges)
                abs_charges = [abs(c) for c in charges]
                result["max_abs_partial_charge"] = max(abs_charges)
                result["min_abs_partial_charge"] = min(abs_charges)
        except Exception as e:
            logger.warning(f"部分電荷の計算に失敗しました: {str(e)}")

        # EState指数
        try:
            estate_indices = EState.EStateIndices(mol)
            if estate_indices:
                result["max_estate_index"] = max(estate_indices)
                result["min_estate_index"] = min(estate_indices)
                abs_estate = [abs(i) for i in estate_indices]
                result["max_abs_estate_index"] = max(abs_estate)
                result["min_abs_estate_index"] = min(abs_estate)
        except Exception as e:
            logger.warning(f"EState指数の計算に失敗しました: {str(e)}")
            
        # グラフ指数
        try:
            result["balaban_j"] = GraphDescriptors.BalabanJ(mol)
            result["bertz_ct"] = GraphDescriptors.BertzCT(mol)
            result["ipc"] = GraphDescriptors.Ipc(mol)
            result["hall_kier_alpha"] = GraphDescriptors.HallKierAlpha(mol)
            result["kappa1"] = GraphDescriptors.Kappa1(mol)
            result["kappa2"] = GraphDescriptors.Kappa2(mol)
            result["kappa3"] = GraphDescriptors.Kappa3(mol)
            result["chi0"] = GraphDescriptors.Chi0(mol)
            result["chi1"] = GraphDescriptors.Chi1(mol)
            result["chi0v"] = GraphDescriptors.Chi0v(mol)
            result["chi1v"] = GraphDescriptors.Chi1v(mol)
        except Exception as e:
            logger.warning(f"グラフ指数の計算に失敗しました: {str(e)}")

        # QED薬剤様性
        try:
            result["qed"] = QED.qed(mol)
        except Exception as e:
            logger.warning(f"QEDの計算に失敗しました: {str(e)}")
        
        # フラグメント解析（官能基カウント）
        for name in dir(Fragments):
            if name.startswith('fr_'):
                try:
                    func = getattr(Fragments, name)
                    if callable(func):
                        result[name] = func(mol)
                except Exception as e:
                    logger.debug(f"{name}の計算に失敗しました: {str(e)}")
        
        # フィルター判定
        
        # Lipinski's Rule of Five
        result["lipinski_molecular_weight_ok"] = result["molecular_weight"] <= 500
        result["lipinski_logp_ok"] = result["logp"] <= 5
        result["lipinski_h_donors_ok"] = result["num_h_donors"] <= 5
        result["lipinski_h_acceptors_ok"] = result["num_h_acceptors"] <= 10
        result["lipinski_pass"] = (
            result["lipinski_molecular_weight_ok"] and
            result["lipinski_logp_ok"] and
            result["lipinski_h_donors_ok"] and
            result["lipinski_h_acceptors_ok"]
        )
        
        # Veber's Rules
        result["veber_rotatable_bonds_ok"] = result["num_rotatable_bonds"] <= 10
        result["veber_tpsa_ok"] = result["tpsa"] <= 140
        result["veber_pass"] = result["veber_rotatable_bonds_ok"] and result["veber_tpsa_ok"]
        
        # Ghose filter
        result["ghose_molecular_weight_ok"] = 160 <= result["molecular_weight"] <= 480
        result["ghose_logp_ok"] = -0.4 <= result["logp"] <= 5.6
        result["ghose_atom_count_ok"] = 20 <= result["heavy_atom_count"] <= 70
        result["ghose_molar_refractivity_ok"] = 40 <= result["mol_mr"] <= 130
        result["ghose_pass"] = (
            result["ghose_molecular_weight_ok"] and
            result["ghose_logp_ok"] and
            result["ghose_atom_count_ok"] and
            result["ghose_molar_refractivity_ok"]
        )
        
        # Egan filter
        result["egan_logp_ok"] = result["logp"] <= 5.88
        result["egan_tpsa_ok"] = result["tpsa"] <= 131.6
        result["egan_pass"] = result["egan_logp_ok"] and result["egan_tpsa_ok"]
        
        # Muegge filter
        result["muegge_molecular_weight_ok"] = 200 <= result["molecular_weight"] <= 600
        result["muegge_logp_ok"] = -2 <= result["logp"] <= 5
        result["muegge_tpsa_ok"] = result["tpsa"] <= 150
        result["muegge_ring_count_ok"] = result["ring_count"] <= 7
        result["muegge_h_acceptors_ok"] = result["num_h_acceptors"] <= 10
        result["muegge_h_donors_ok"] = result["num_h_donors"] <= 5
        result["muegge_rotatable_bonds_ok"] = result["num_rotatable_bonds"] < 15
        result["muegge_pass"] = (
            result["muegge_molecular_weight_ok"] and
            result["muegge_logp_ok"] and
            result["muegge_tpsa_ok"] and
            result["muegge_ring_count_ok"] and
            result["muegge_h_acceptors_ok"] and
            result["muegge_h_donors_ok"] and
            result["muegge_rotatable_bonds_ok"]
        )
        
        # PAINS filter
        result["pains_free"] = True
        result["pains_num_alerts"] = 0
        result["pains_alerts"] = []
        
        try:
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
                result["pains_num_alerts"] = len(matches)
                
                # マッチしたアラートの詳細情報を取得
                pains_alerts = []
                for match in matches:
                    pains_alerts.append({
                        "description": match.GetDescription(),
                        "smarts": match.GetSmarts() if hasattr(match, "GetSmarts") else None
                    })
                result["pains_alerts"] = pains_alerts
        except Exception as e:
            logger.warning(f"PAINSフィルターの適用中にエラーが発生しました: {str(e)}")
        
        # 総合判定
        result["all_filters_passed"] = (
            result["lipinski_pass"] and
            result["veber_pass"] and
            result["ghose_pass"] and
            result["egan_pass"] and
            result["muegge_pass"] and
            result["pains_free"]
        )
        
    except Exception as e:
        logger.error(f"分子特性の計算中にエラーが発生しました: {str(e)}")
    
    return result