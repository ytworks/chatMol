"""
chatMol - Library for calculating molecular properties from molecular structures
"""

from .properties import calculate_molecular_features, get_property_descriptions, get_available_properties, get_feature_descriptions
from .io import add_properties_to_dataframe

__version__ = "0.1.0"