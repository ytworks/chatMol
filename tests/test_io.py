"""
Test module for the chatMol io module.
Tests the input/output functionality for molecular data.
"""
import pytest
import pandas as pd
from chatmol.io import add_properties_to_dataframe
from chatmol.properties import calculate_molecular_features

# テストデータ
TEST_MOLECULES = [
    {"compound_id": "ASP", "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"},  # アスピリン
    {"compound_id": "PCM", "smiles": "CC(=O)NC1=CC=C(C=C1)O"},  # パラセタモール
    {"compound_id": "IBP", "smiles": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"}  # イブプロフェン
]


class TestIOFunctionality:
    """Test class for input/output functionality."""
    
    def test_add_properties_to_dataframe(self):
        """Test adding molecular properties to a DataFrame."""
        # テスト用のデータフレームを作成
        df = pd.DataFrame(TEST_MOLECULES)
        
        # 分子特性を計算
        feature_results = [
            calculate_molecular_features(mol["smiles"]) for mol in TEST_MOLECULES
        ]
        
        # データフレームに特性を追加
        add_properties_to_dataframe(df, feature_results)
        
        # 主要なプロパティがデータフレームに追加されていることを確認
        assert "molecular_weight" in df.columns
        assert "logp" in df.columns
        assert "formula" in df.columns
        assert "num_h_donors" in df.columns
        assert "num_h_acceptors" in df.columns
        assert "tpsa" in df.columns
        
        # データの型と値が期待通りであることを確認（アスピリンの例）
        assert isinstance(df.loc[0, "molecular_weight"], float)
        assert round(df.loc[0, "molecular_weight"], 1) == 180.2  # アスピリンの分子量
        assert isinstance(df.loc[0, "formula"], str)
        assert df.loc[0, "formula"] == "C9H8O4"  # アスピリンの分子式
    
    def test_column_name_conflict_resolution(self):
        """Test handling of column name conflicts when adding properties."""
        # 既に'molecular_weight'カラムがあるデータフレームを作成
        df = pd.DataFrame({
            "compound_id": ["ASP", "PCM", "IBP"],
            "smiles": [mol["smiles"] for mol in TEST_MOLECULES],
            "molecular_weight": [100, 200, 300]  # ダミーの値
        })
        
        # 分子特性を計算
        feature_results = [
            calculate_molecular_features(mol["smiles"]) for mol in TEST_MOLECULES
        ]
        
        # データフレームに特性を追加
        add_properties_to_dataframe(df, feature_results)
        
        # カラム名の競合が適切に解決されていることを確認
        assert "molecular_weight" in df.columns  # 元のカラムは保持
        assert "molecular_weight_calculated" in df.columns  # 新しいカラムは名前が変更されている
        
        # 元の値とRDKitで計算された値が異なることを確認
        assert df.loc[0, "molecular_weight"] == 100  # 元の値
        assert round(df.loc[0, "molecular_weight_calculated"], 1) == 180.2  # 計算された値（アスピリン）
    
    def test_empty_dataframe(self):
        """Test handling of empty DataFrames."""
        # 空のデータフレームを作成
        df = pd.DataFrame(columns=["compound_id", "smiles"])
        
        # 空の結果リスト
        feature_results = []
        
        # データフレームに特性を追加
        add_properties_to_dataframe(df, feature_results)
        
        # 元のカラムは保持されていることを確認
        assert "compound_id" in df.columns
        assert "smiles" in df.columns
        assert len(df) == 0  # 行数は0のまま
    
    def test_partial_property_results(self):
        """Test handling of partial property results (some molecules have properties others don't)."""
        # テスト用のデータフレームを作成
        df = pd.DataFrame({
            "compound_id": ["ASP", "INVALID", "IBP"],
            "smiles": [
                "CC(=O)OC1=CC=CC=C1C(=O)O",  # アスピリン
                "invalid_smiles",             # 無効なSMILES
                "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"  # イブプロフェン
            ]
        })
        
        # 分子特性を計算（2番目は無効なSMILES）
        feature_results = [
            calculate_molecular_features(df.loc[0, "smiles"]),
            calculate_molecular_features(df.loc[1, "smiles"]),
            calculate_molecular_features(df.loc[2, "smiles"])
        ]
        
        # データフレームに特性を追加
        add_properties_to_dataframe(df, feature_results)
        
        # 正常なSMILESに対しては値が入っていることを確認
        assert "molecular_weight" in df.columns
        assert pd.notna(df.loc[0, "molecular_weight"])
        assert pd.notna(df.loc[2, "molecular_weight"])
        
        # 無効なSMILESに対してはNoneまたはNaNになっていることを確認
        assert pd.isna(df.loc[1, "molecular_weight"])