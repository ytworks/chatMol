"""
chatMol - 分子構造からの物性値計算を行うライブラリ
"""

from .properties import calculate_molecular_weight, calculate_properties
from .io import process_csv_data

__version__ = "0.1.0"