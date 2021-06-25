import bz2
import gzip
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import _pickle as cPickle
import io
import sys

from collections import Counter
from PIL import Image

from rdkit import Chem, RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import AllChem, Draw, PandasTools, rdFMCS, rdMolDescriptors
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.Draw import IPythonConsole, MolsToGridImage
from rdkit.DataStructs import cDataStructs
from rdkit_ipynb_tools import tools, pipeline as p


# Set file path and format
file = 'data/complete_chemical_data.sdf.gz'

# Read SDF
moldf = PandasTools.LoadSDF(file, molColName='Mol', smilesName='SMILES');
print('Original data: ', moldf.shape)

# Remove missing RDKit molecules
moldf = moldf[pd.notnull(moldf['Mol'])]
if 'StandardizerResult' in moldf.columns:
    moldf = moldf.drop(columns='StandardizerResult')

moldf.replace(r'^\s*$', np.nan, regex=True, inplace=True)
moldf.drop(columns='ID', inplace=True)
# Columns
print('Kept data: ', moldf.shape)
moldf.head(1)
from molvs.validate import Validator
fmt = '%(asctime)s - %(levelname)s - %(validation)s - %(message)s'
validator = Validator(log_format=fmt)
print('\n Problematic structures: \n', validator.validate(moldf))