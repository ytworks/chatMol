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
        Dict: Dictionary containing property names as keys and their description information
    """
    return {
        # Basic properties
        "molecular_weight": {
            "ja": "分子量",
            "en": "Molecular Weight",
            "description_ja": "分子の平均分子量（各元素の平均原子量に基づく）。標準的な分子量を表す。",
            "description_en": "Average molecular weight of the molecule based on the average atomic masses of elements. Represents the standard molecular weight.",
            "module": "Descriptors.MolWt"
        },
        "exact_mol_wt": {
            "ja": "厳密分子量",
            "en": "Exact Molecular Weight",
            "description_ja": "同位体組成を考慮した厳密な分子量（モノアイソトピック質量）。小数点まで厳密に算出された質量。",
            "description_en": "Exact molecular weight considering isotopic composition (monoisotopic mass). Mass calculated precisely to decimal places.",
            "module": "Descriptors.ExactMolWt"
        },
        "heavy_atom_mol_wt": {
            "ja": "重原子分子量",
            "en": "Heavy Atom Molecular Weight",
            "description_ja": "水素を無視した分子の平均分子量。炭素など重原子のみの質量合計。",
            "description_en": "Average molecular weight ignoring hydrogens. Sum of the mass of heavy atoms only, such as carbon.",
            "module": "Descriptors.HeavyAtomMolWt"
        },
        "formula": {
            "ja": "分子式",
            "en": "Molecular Formula",
            "description_ja": "分子の化学式。分子を構成する原子の種類と数を表す。",
            "description_en": "Chemical formula of the molecule. Represents the types and numbers of atoms constituting the molecule.",
            "module": "rdMolDescriptors.CalcMolFormula"
        },
        
        # Lipophilicity/Hydrophilicity
        "logp": {
            "ja": "オクタノール/水分配係数",
            "en": "LogP",
            "description_ja": "分配係数LogP（1-オクタノールと水の間の分配係数の対数）。Wildman–Crippen法による推算値。疎水性の指標。",
            "description_en": "Partition coefficient LogP (logarithm of the partition coefficient between 1-octanol and water). Estimated by Wildman-Crippen method. An indicator of hydrophobicity.",
            "module": "Descriptors.MolLogP"
        },
        "mol_mr": {
            "ja": "モル屈折率",
            "en": "Molar Refractivity",
            "description_ja": "モル屈折率。Wildman–Crippen法による推算値で、分子の分極率に関連する指標。",
            "description_en": "Molar refractivity. Estimated by Wildman-Crippen method, related to the polarizability of the molecule.",
            "module": "Descriptors.MolMR"
        },
        "tpsa": {
            "ja": "トポロジカル極性表面積",
            "en": "Topological Polar Surface Area",
            "description_ja": "分子中の極性原子（主にOとN、およびそれに結合したH）の表面積の合計。極性表面積が大きいほど細胞膜透過性は低下する傾向がある。",
            "description_en": "Sum of the surface areas of polar atoms in a molecule (mainly O, N, and their bonded H atoms). Larger polar surface area tends to decrease cell membrane permeability.",
            "module": "Descriptors.TPSA"
        },

        # Surface areas
        "labute_asa": {
            "ja": "Labute推定表面積",
            "en": "Labute's Approx. Surface Area",
            "description_ja": "Labuteによる近似分子表面積。MOEソフトウェア由来のアルゴリズムで計算される分子表面積。",
            "description_en": "Approximate molecular surface area by Labute. Molecular surface area calculated using an algorithm derived from MOE software.",
            "module": "Descriptors.LabuteASA"
        },
        
        # H-bonds and atom counts
        "num_h_donors": {
            "ja": "水素結合ドナー数",
            "en": "#H-Bond Donors",
            "description_ja": "分子内の水素結合供与体の個数。一般に–OHや–NH基の数（Lipinskiの定義に基づく）。",
            "description_en": "Number of hydrogen bond donors in the molecule. Generally the count of -OH and -NH groups (based on Lipinski's definition).",
            "module": "Descriptors.NumHDonors"
        },
        "num_h_acceptors": {
            "ja": "水素結合アクセプター数",
            "en": "#H-Bond Acceptors",
            "description_ja": "分子内の水素結合受容体の個数。一般にOやN原子の数（Lipinskiの定義に基づく）。",
            "description_en": "Number of hydrogen bond acceptors in the molecule. Generally the count of O and N atoms (based on Lipinski's definition).",
            "module": "Descriptors.NumHAcceptors"
        },
        "num_rotatable_bonds": {
            "ja": "回転可能結合数",
            "en": "#Rotatable Bonds",
            "description_ja": "分子内の回転可能な単結合の数（末端の単結合や二重結合環境を除く）。分子の柔軟性の指標。",
            "description_en": "Number of rotatable single bonds in the molecule (excluding terminal single bonds and those in double bond environments). An indicator of molecular flexibility.",
            "module": "Descriptors.NumRotatableBonds"
        },
        "heavy_atom_count": {
            "ja": "重原子数",
            "en": "#Heavy Atoms",
            "description_ja": "分子中の重原子（非水素原子）の数。分子サイズの指標の一つ。",
            "description_en": "Number of heavy atoms (non-hydrogen atoms) in the molecule. One of the indicators of molecular size.",
            "module": "Descriptors.HeavyAtomCount"
        },
        "num_hetero_atoms": {
            "ja": "ヘテロ原子数",
            "en": "#Heteroatoms",
            "description_ja": "分子中のヘテロ原子（炭素以外の元素）の数。例えばN, O, Sなどの個数。",
            "description_en": "Number of heteroatoms (elements other than carbon) in the molecule. For example, count of N, O, S atoms.",
            "module": "Descriptors.NumHeteroatoms"
        },
        "no_count": {
            "ja": "N/O原子数",
            "en": "#Nitrogen/Oxygen Atoms",
            "description_ja": "分子中の窒素原子および酸素原子の合計数。LipinskiのH受容体数と近似的に対応。",
            "description_en": "Total number of nitrogen and oxygen atoms in the molecule. Approximately corresponds to Lipinski's H-acceptor count.",
            "module": "Lipinski.NOCount"
        },
        "nhoh_count": {
            "ja": "OH/NH基数",
            "en": "#OH/NH Groups",
            "description_ja": "分子中のヒドロキシ基(-OH)および一次・二次アミン基(-NH)の数。LipinskiのHドナー数と近い概念。",
            "description_en": "Number of hydroxyl (-OH) and primary/secondary amine (-NH) groups in the molecule. Similar concept to Lipinski's H-donor count.",
            "module": "Lipinski.NHOHCount"
        },
        "num_valence_electrons": {
            "ja": "原子価電子数",
            "en": "#Valence Electrons",
            "description_ja": "分子内の全原子の価電子の総数。分子全体の価電子の数で、電荷や元素組成の指標。",
            "description_en": "Total number of valence electrons of all atoms in the molecule. Count of valence electrons for the entire molecule, an indicator of charge and elemental composition.",
            "module": "Descriptors.NumValenceElectrons"
        },
        
        # Ring information
        "num_aromatic_rings": {
            "ja": "芳香環数",
            "en": "#Aromatic Rings",
            "description_ja": "分子内の芳香族環の数。ベンゼン環など芳香環構造の個数。",
            "description_en": "Number of aromatic rings in the molecule. Count of aromatic ring structures such as benzene rings.",
            "module": "Descriptors.NumAromaticRings"
        },
        "num_aliphatic_rings": {
            "ja": "脂肪族環数",
            "en": "#Aliphatic Rings",
            "description_ja": "分子内の脂肪族環の数。非芳香族の環構造（シクロアルカンなど）の個数。",
            "description_en": "Number of aliphatic rings in the molecule. Count of non-aromatic ring structures (such as cycloalkanes).",
            "module": "Descriptors.NumAliphaticRings"
        },
        "num_saturated_rings": {
            "ja": "飽和環数",
            "en": "#Saturated Rings",
            "description_ja": "分子内の飽和環の数。二重結合を含まない環（飽和脂肪族環）の個数。",
            "description_en": "Number of saturated rings in the molecule. Count of rings without double bonds (saturated aliphatic rings).",
            "module": "Descriptors.NumSaturatedRings"
        },
        "num_aromatic_carbocycles": {
            "ja": "芳香族炭素環数",
            "en": "#Aromatic Carbocycles",
            "description_ja": "芳香族炭素環（全原子が炭素の芳香環）の数。",
            "description_en": "Number of aromatic carbocycles (aromatic rings where all atoms are carbon).",
            "module": "Descriptors.NumAromaticCarbocycles"
        },
        "num_aromatic_heterocycles": {
            "ja": "芳香族複素環数",
            "en": "#Aromatic Heterocycles",
            "description_ja": "芳香族複素環（炭素以外の原子を含む芳香環）の数。",
            "description_en": "Number of aromatic heterocycles (aromatic rings containing atoms other than carbon).",
            "module": "Descriptors.NumAromaticHeterocycles"
        },
        "num_aliphatic_carbocycles": {
            "ja": "脂肪族炭素環数",
            "en": "#Aliphatic Carbocycles",
            "description_ja": "脂肪族炭素環（炭素のみからなる非芳香環）の数。",
            "description_en": "Number of aliphatic carbocycles (non-aromatic rings consisting only of carbon atoms).",
            "module": "Descriptors.NumAliphaticCarbocycles"
        },
        "num_aliphatic_heterocycles": {
            "ja": "脂肪族複素環数",
            "en": "#Aliphatic Heterocycles",
            "description_ja": "脂肪族複素環（ヘテロ原子を含む非芳香環）の数。",
            "description_en": "Number of aliphatic heterocycles (non-aromatic rings containing heteroatoms).",
            "module": "Descriptors.NumAliphaticHeterocycles"
        },
        "num_saturated_carbocycles": {
            "ja": "飽和炭素環数",
            "en": "#Saturated Carbocycles",
            "description_ja": "飽和炭素環（完全に単結合のみからなる炭素環）の数。",
            "description_en": "Number of saturated carbocycles (carbon rings consisting entirely of single bonds).",
            "module": "Descriptors.NumSaturatedCarbocycles"
        },
        "num_saturated_heterocycles": {
            "ja": "飽和複素環数",
            "en": "#Saturated Heterocycles",
            "description_ja": "飽和複素環（完全に単結合のみからなる複素環）の数。",
            "description_en": "Number of saturated heterocycles (heterocyclic rings consisting entirely of single bonds).",
            "module": "Descriptors.NumSaturatedHeterocycles"
        },
        "ring_count": {
            "ja": "環数（総数）",
            "en": "Ring Count",
            "description_ja": "分子内の環構造の総数。全ての大小の環をカウントしたもの。",
            "description_en": "Total number of ring structures in the molecule. Count of all rings regardless of size.",
            "module": "Descriptors.RingCount"
        },
        
        # Bond/functional group counts
        "num_amide_bonds": {
            "ja": "アミド結合数",
            "en": "#Amide Bonds",
            "description_ja": "分子中のアミド結合（–C(=O)N–結合）の数。",
            "description_en": "Number of amide bonds (-C(=O)N- bonds) in the molecule.",
            "module": "Descriptors.NumAmideBonds"
        },
        "fraction_csp3": {
            "ja": "炭素sp³割合",
            "en": "Fraction of C sp³",
            "description_ja": "炭素原子のうちsp³混成状態にあるものの割合。値域0〜1で、分子の飽和度を示す指標。",
            "description_en": "Fraction of carbon atoms in sp³ hybridization state. Range 0-1, an indicator of molecular saturation degree.",
            "module": "Descriptors.FractionCSP3"
        },
        
        # Molecular complexity
        "num_spiro_atoms": {
            "ja": "スピロ原子数",
            "en": "#Spiro Atoms",
            "description_ja": "スピロ原子の数。2つの環が1つの原子を共有して融合した構造（スピロ環）における共有原子の個数。",
            "description_en": "Number of spiro atoms. Count of shared atoms in structures where two rings are fused by sharing a single atom (spiro rings).",
            "module": "Descriptors.NumSpiroAtoms"
        },
        "num_bridgehead_atoms": {
            "ja": "ブリッジヘッド原子数",
            "en": "#Bridgehead Atoms",
            "description_ja": "ブリッジヘッド原子の数。複数の環が少なくとも2つの結合を共有して融合した構造（ブリッジヘッド）における共有原子の個数。",
            "description_en": "Number of bridgehead atoms. Count of shared atoms in structures where multiple rings are fused by sharing at least two bonds (bridgehead structures).",
            "module": "Descriptors.NumBridgeheadAtoms"
        },
        "num_stereo_centers": {
            "ja": "立体中心数",
            "en": "#Stereo Centers",
            "description_ja": "分子内の立体中心（キラル中心）となっている原子の数。明確に(R/S等)指定された立体中心の個数。",
            "description_en": "Number of stereogenic centers (chiral centers) in the molecule. Count of atoms with clearly specified (R/S, etc.) stereochemistry.",
            "module": "Descriptors.NumAtomStereoCenters"
        },
        "num_unspecified_stereo_centers": {
            "ja": "未指定立体中心数",
            "en": "#Unspecified Stereo Centers",
            "description_ja": "分子内の立体中心のうち立体配置が未定義のものの数（立体化学が指定されていないキラル原子の個数）。",
            "description_en": "Number of stereogenic centers with undefined stereochemistry (count of chiral atoms without specified stereochemistry).",
            "module": "Descriptors.NumUnspecifiedAtomStereoCenters"
        },
        
        # Graph indices
        "balaban_j": {
            "ja": "BalabanのJ指数",
            "en": "Balaban's J Index",
            "description_ja": "Balabanが提唱した分子接続指数J。分子の接続グラフから計算されるトップロジカル指数で、分子の密結合度合いを表す。",
            "description_en": "Balaban's molecular connectivity index J. A topological index calculated from the molecular connection graph, representing the degree of molecular connectivity.",
            "module": "GraphDescriptors.BalabanJ"
        },
        "bertz_ct": {
            "ja": "Bertzの複雑度指数",
            "en": "Bertz's Complexity Index",
            "description_ja": "Bertzによる分子構造の複雑さ（Complexity）を定量化する指数。原子と結合の組み合わせの情報量に基づき分子の複雑度を評価する。",
            "description_en": "Bertz's index quantifying molecular structure complexity. Evaluates molecular complexity based on the information content of atom and bond combinations.",
            "module": "GraphDescriptors.BertzCT"
        },
        "ipc": {
            "ja": "Ipc情報指数",
            "en": "Ipc (Information Content)",
            "description_ja": "分子グラフの情報量記述子Ipc。隣接行列の特性多項式係数から算出される情報量に基づく指標。",
            "description_en": "Information content descriptor Ipc of the molecular graph. An index based on information calculated from characteristic polynomial coefficients of the adjacency matrix.",
            "module": "GraphDescriptors.Ipc"
        },
        "hall_kier_alpha": {
            "ja": "Hall–Kierのαパラメータ",
            "en": "Hall-Kier Alpha",
            "description_ja": "Kierらによる分子の補正パラメータα。一般に分子のサイズや分枝に関する補正項で、κ形状指数の計算にも用いられる。",
            "description_en": "Correction parameter α for molecules by Kier et al. Generally a correction term related to molecular size and branching, also used in calculating κ shape indices.",
            "module": "GraphDescriptors.HallKierAlpha"
        },

        # Kappa shape indices
        "kappa1": {
            "ja": "κ形状指数1",
            "en": "Kappa Shape Index 1",
            "description_ja": "Kierの形状指数。与えられた原子数の中で最も線状な場合および最も分岐した場合を基準に、分子の形状を表す指数。Kappa1は分子の分岐度を表す。",
            "description_en": "Kier's shape index. An index representing molecular shape based on the most linear and most branched cases for a given number of atoms. Kappa1 represents the degree of molecular branching.",
            "module": "GraphDescriptors.Kappa1"
        },
        "kappa2": {
            "ja": "κ形状指数2",
            "en": "Kappa Shape Index 2",
            "description_ja": "Kierの形状指数。Kappa2は分子の空間的な広がり（面状）を表す。",
            "description_en": "Kier's shape index. Kappa2 represents the spatial extent (planar aspect) of the molecule.",
            "module": "GraphDescriptors.Kappa2"
        },
        "kappa3": {
            "ja": "κ形状指数3",
            "en": "Kappa Shape Index 3",
            "description_ja": "Kierの形状指数。Kappa3は分子の空間的な広がり（立体的）を表す。",
            "description_en": "Kier's shape index. Kappa3 represents the spatial extent (three-dimensional aspect) of the molecule.",
            "module": "GraphDescriptors.Kappa3"
        },
        
        # Chi connectivity indices
        "chi0": {
            "ja": "χ0接続性指数",
            "en": "Chi0 Connectivity Index",
            "description_ja": "KierとHallが提案した分子接続性指数。Chi0は0次の接続性指数。",
            "description_en": "Molecular connectivity index proposed by Kier and Hall. Chi0 is the zeroth-order connectivity index.",
            "module": "GraphDescriptors.Chi0"
        },
        "chi1": {
            "ja": "χ1接続性指数",
            "en": "Chi1 Connectivity Index",
            "description_ja": "KierとHallが提案した分子接続性指数。Chi1は1次の接続性指数。",
            "description_en": "Molecular connectivity index proposed by Kier and Hall. Chi1 is the first-order connectivity index.",
            "module": "GraphDescriptors.Chi1"
        },
        "chi0v": {
            "ja": "χ0v接続性指数（原子価考慮）",
            "en": "Chi0v Connectivity Index",
            "description_ja": "KierとHallが提案した分子接続性指数。Chi0vは原子価を考慮した0次の接続性指数。",
            "description_en": "Molecular connectivity index proposed by Kier and Hall. Chi0v is the zeroth-order connectivity index considering valence.",
            "module": "GraphDescriptors.Chi0v"
        },
        "chi1v": {
            "ja": "χ1v接続性指数（原子価考慮）",
            "en": "Chi1v Connectivity Index",
            "description_ja": "KierとHallが提案した分子接続性指数。Chi1vは原子価を考慮した1次の接続性指数。",
            "description_en": "Molecular connectivity index proposed by Kier and Hall. Chi1v is the first-order connectivity index considering valence.",
            "module": "GraphDescriptors.Chi1v"
        },
        
        # QED drug-likeness
        "qed": {
            "ja": "QED薬剤様性スコア",
            "en": "QED Drug-Likeness Score",
            "description_ja": "QED (Quantitative Estimation of Drug-likeness)。分子量、LogP、TPSA、Hドナー/アクセプター数、芳香環数、回転結合数など複数の性質を総合して0〜1の範囲で算出される薬剤様性指標。値が高いほどドラッグライクとされる。",
            "description_en": "QED (Quantitative Estimation of Drug-likeness). A drug-likeness indicator calculated in the range of 0-1 by combining multiple properties such as molecular weight, LogP, TPSA, H-donor/acceptor counts, aromatic ring count, and rotatable bond count. Higher values indicate more drug-like properties.",
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
    Get descriptions of all molecular features (properties and filters)
    
    Returns:
        Dict: Dictionary containing feature descriptions
    """
    # Get property descriptions
    feature_descriptions = get_property_descriptions()
    
    # Add filter information in the same format
    for filter_name, filter_info in MOLECULAR_FILTERS.items():
        # Convert to the same format as property descriptions
        feature_descriptions[filter_name] = {
            "ja": filter_info.get("ja", filter_name),  # Use key name if Japanese name not available
            "en": filter_info.get("en", filter_name),  # Use key name if English name not available
            "description_ja": filter_info.get("description_ja", ""),
            "description_en": filter_info.get("description_en", filter_info.get("description", "")),
            "is_filter": True,  # Flag to indicate this is a filter
            "result_key": filter_info.get("result_key", "")  # Filter-specific information
        }
    
    return feature_descriptions


# Filter definitions with their result keys and descriptions in both English and Japanese
MOLECULAR_FILTERS = {
    "lipinski": {
        "result_key": "all_rules_passed",
        "ja": "リピンスキーの5原則",
        "en": "Lipinski's Rule of Five",
        "description_en": "Lipinski's Rule of Five (MW≤500, LogP≤5, HBD≤5, HBA≤10). Rules for drug-like compounds with good oral bioavailability.",
        "description_ja": "リピンスキーの5原則（分子量≤500、LogP≤5、水素結合ドナー≤5、水素結合アクセプター≤10）。経口投与可能な薬剤様化合物のための指標。"
    },
    "veber": {
        "result_key": "all_rules_passed",
        "ja": "ヴェーバーの法則",
        "en": "Veber's Rules",
        "description_en": "Veber's rules (TPSA≤140 Å², RotBonds≤10). Additional criteria for good oral bioavailability in rats.",
        "description_ja": "ヴェーバーの法則（極性表面積≤140 Å²、回転可能結合≤10）。ラットにおける良好な経口バイオアベイラビリティのための追加基準。"
    },
    "ghose": {
        "result_key": "all_rules_passed",
        "ja": "ゴーシュフィルター",
        "en": "Ghose Filter",
        "description_en": "Ghose filter (160≤MW≤480, -0.4≤LogP≤5.6, 20≤atoms≤70, 40≤MR≤130). Physicochemical parameters for drug-likeness.",
        "description_ja": "ゴーシュフィルター（160≤分子量≤480、-0.4≤LogP≤5.6、20≤原子数≤70、40≤モル屈折率≤130）。薬剤様性のための物理化学的パラメータ。"
    },
    "egan": {
        "result_key": "all_rules_passed",
        "ja": "イーガンフィルター",
        "en": "Egan Filter",
        "description_en": "Egan filter (LogP≤5.88, TPSA≤131.6). Absorption prediction based on passive intestinal absorption.",
        "description_ja": "イーガンフィルター（LogP≤5.88、極性表面積≤131.6）。受動的な腸管吸収に基づく吸収性予測。"
    },
    "muegge": {
        "result_key": "all_rules_passed",
        "ja": "ミュッゲフィルター",
        "en": "Muegge Filter",
        "description_en": "Muegge filter (200≤MW≤600, -2≤LogP≤5, TPSA≤150, rings≤7, HBA≤10, HBD≤5, RotBonds<15). Pharmacophore points for good oral bioavailability.",
        "description_ja": "ミュッゲフィルター（200≤分子量≤600、-2≤LogP≤5、極性表面積≤150、環数≤7、水素結合アクセプター≤10、水素結合ドナー≤5、回転可能結合<15）。良好な経口バイオアベイラビリティのためのファーマコフォアポイント。"
    },
    "pains": {
        "result_key": "pains_free",
        "ja": "PAINSフィルター",
        "en": "PAINS Filter",
        "description_en": "PAINS filter (screens for pan-assay interference compounds). Identifies compounds likely to cause false positives in biological assays.",
        "description_ja": "PAINSフィルター（パンアッセイ干渉化合物のスクリーニング）。生物学的アッセイで偽陽性を引き起こす可能性の高い化合物を特定する。"
    }
}


def calculate_molecular_features(
    smiles: str, 
    use_rdkit_mol: bool = False
) -> Dict[str, Any]:
    """
    Calculates all molecular properties and filter evaluations for a molecule and returns them in a flat dictionary
    
    Args:
        smiles: Molecular structure in SMILES notation
        use_rdkit_mol: If True, also returns the RDKit molecule object (for reuse)
    
    Returns:
        Dict: Dictionary containing all calculated properties and filter results in a flat structure
    """
    # Initialize result dictionary (flat structure)
    result = {
        "smiles": smiles
    }
    
    # Check if RDKit is available
    if not rdkit_available:
        logger.error("RDKit is not installed. Cannot calculate molecular properties.")
        if use_rdkit_mol:
            result["mol"] = None
        return result
    
    try:
        # Create RDKit molecule object from SMILES string
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Invalid SMILES string: {smiles}")
            if use_rdkit_mol:
                result["mol"] = None
            return result
        
        # Save molecule object (optional)
        if use_rdkit_mol:
            result["mol"] = mol
        
        # Basic molecular properties
        result["molecular_weight"] = Descriptors.MolWt(mol)
        result["exact_mol_wt"] = Descriptors.ExactMolWt(mol)
        result["heavy_atom_mol_wt"] = Descriptors.HeavyAtomMolWt(mol)
        result["formula"] = rdMolDescriptors.CalcMolFormula(mol)
        
        # Lipophilicity/Hydrophilicity
        result["logp"] = Descriptors.MolLogP(mol)
        result["mol_mr"] = Descriptors.MolMR(mol)
        result["tpsa"] = Descriptors.TPSA(mol)
        
        # Surface area
        result["labute_asa"] = Descriptors.LabuteASA(mol)
        
        # H-bonds and atom counts
        result["num_h_donors"] = Descriptors.NumHDonors(mol)
        result["num_h_acceptors"] = Descriptors.NumHAcceptors(mol)
        result["num_rotatable_bonds"] = Descriptors.NumRotatableBonds(mol)
        result["heavy_atom_count"] = Descriptors.HeavyAtomCount(mol)
        result["num_hetero_atoms"] = Descriptors.NumHeteroatoms(mol)
        result["no_count"] = Lipinski.NOCount(mol)
        result["nhoh_count"] = Lipinski.NHOHCount(mol)
        result["num_valence_electrons"] = Descriptors.NumValenceElectrons(mol)
        
        # Ring information
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
        
        # Bond/functional group counts
        result["num_amide_bonds"] = Descriptors.NumAmideBonds(mol)
        result["fraction_csp3"] = Descriptors.FractionCSP3(mol)
        
        # Molecular complexity
        result["num_spiro_atoms"] = Descriptors.NumSpiroAtoms(mol)
        result["num_bridgehead_atoms"] = Descriptors.NumBridgeheadAtoms(mol)
        result["num_stereo_centers"] = Descriptors.NumAtomStereoCenters(mol)
        result["num_unspecified_stereo_centers"] = Descriptors.NumUnspecifiedAtomStereoCenters(mol)
        
        # Charge-related properties
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
            logger.warning(f"Failed to calculate partial charges: {str(e)}")

        # EState indices
        try:
            estate_indices = EState.EStateIndices(mol)
            if estate_indices:
                result["max_estate_index"] = max(estate_indices)
                result["min_estate_index"] = min(estate_indices)
                abs_estate = [abs(i) for i in estate_indices]
                result["max_abs_estate_index"] = max(abs_estate)
                result["min_abs_estate_index"] = min(abs_estate)
        except Exception as e:
            logger.warning(f"Failed to calculate EState indices: {str(e)}")
            
        # Graph indices
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
            logger.warning(f"Failed to calculate graph indices: {str(e)}")

        # QED drug-likeness
        try:
            result["qed"] = QED.qed(mol)
        except Exception as e:
            logger.warning(f"Failed to calculate QED: {str(e)}")
        
        # Fragment analysis (functional group counts)
        for name in dir(Fragments):
            if name.startswith('fr_'):
                try:
                    func = getattr(Fragments, name)
                    if callable(func):
                        result[name] = func(mol)
                except Exception as e:
                    logger.debug(f"Failed to calculate {name}: {str(e)}")
        
        # Filter evaluations
        
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
            # Create filter catalog for PAINS
            params = FilterCatalogParams()
            params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
            params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
            params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
            catalog = FilterCatalog(params)
            
            # Check if the molecule has PAINS patterns
            if catalog.HasMatch(mol):
                # Get matched entries
                matches = catalog.GetMatches(mol)
                result["pains_free"] = False
                result["pains_num_alerts"] = len(matches)
                
                # Get detailed information for matched alerts
                pains_alerts = []
                for match in matches:
                    pains_alerts.append({
                        "description": match.GetDescription(),
                        "smarts": match.GetSmarts() if hasattr(match, "GetSmarts") else None
                    })
                result["pains_alerts"] = pains_alerts
        except Exception as e:
            logger.warning(f"Error occurred while applying PAINS filter: {str(e)}")
        
        # Overall evaluation
        result["all_filters_passed"] = (
            result["lipinski_pass"] and
            result["veber_pass"] and
            result["ghose_pass"] and
            result["egan_pass"] and
            result["muegge_pass"] and
            result["pains_free"]
        )
        
    except Exception as e:
        logger.error(f"Error occurred while calculating molecular properties: {str(e)}")
    
    return result