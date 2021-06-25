import csv
import numpy as np
import pandas as pd
import pickle
import joblib

from cheminformatics import cheminformatics

from chembl_structure_pipeline import standardizer
from molvs.metal import MetalDisconnector

from rdkit import Chem
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import AllChem, PandasTools

from sklearn.metrics.pairwise import pairwise_distances
from sklearn.preprocessing import StandardScaler

from IPython.core.display import display, HTML
display(HTML("<style>.container { width:90% !important; }</style>"))

def warn(*args, **kwargs):
    pass
import warnings
warnings.filterwarnings("ignore")
warnings.warn = warn

# Read training set SDF
file = 'data/complete_chemical_data_molvs.sdf.gz'
moldf = PandasTools.LoadSDF(file, molColName='Mol', smilesName='SMILES')

# Morgan descriptors
morgan_train = pd.read_csv('data/descriptors/morgan.csv')

# MACCS descriptors
maccs_train = pd.read_csv('data/descriptors/maccs.csv')

# AtomPair
atom_pair_train = pd.read_csv('data/descriptors/atom_pair.csv')

# Avalon
avalon_train = pd.read_csv('data/descriptors/avalon.csv')

# Biological data
# DrugMatrix
drug_matrix = pd.read_csv('data/drug_matrix.csv', na_values='?', quoting=csv.QUOTE_NONNUMERIC)

# tg_gates
tg_gates = pd.read_csv('data/tg_gates.csv', na_values='?', quoting=csv.QUOTE_NONNUMERIC)

# list of parameters
with open('param_list.txt', 'r') as f:
    param_list = f.read().split(',')

with open('weight_list.txt', 'r') as f:
    weight_list = f.read().split(',')
    
with open('clin_chem_list.txt', 'r') as f:
    clin_chem_list = f.read().split(',')
    
with open('hemato_list.txt', 'r') as f:
    hemato_list = f.read().split(',')

# Calculate chemical descriptors
def calculate_descriptors(mol):
    # Calculate Morgan
    morgan = cheminformatics.calc_morgan(mol, radius=(2), nBits=(2048))
    morgan = morgan.values.reshape(1,-1)
    
    # Calculate MACCS
    maccs = cheminformatics.calc_maccs(mol)
    maccs = maccs.values.reshape(1,-1)
    
    # Calculate Morgan
    atom_pair = cheminformatics.calc_atom_pair(mol)
    atom_pair = atom_pair.values.reshape(1,-1)
    
    # Calculate Avalon descriptors
    avalon = cheminformatics.calc_avalon(mol, nBits=(2048))
    avalon = avalon.values.reshape(1,-1)
    
    return morgan, maccs, atom_pair, avalon

# Functions to get nearest neighbors in each space
    
def similarity_by_chunk(start, end, X_test, X_train, matrix_len):
    if end > matrix_len:
        end = matrix_len
    return pairwise_distances(X_test, X_train, metric='jaccard', n_jobs=-1)
    
def get_nearest_neighbor(descriptor_train, descriptor_query):
    try:
        X_train = descriptor_train.values
    except:
        X_train = descriptor_train
    X_test = descriptor_query
    matrix_len = X_test.shape[0]
    chunk_size = 5000
    
    nearest_neighbor = []
    val_indexes = []
    dm = []
    similarity = []
    train_indexes = []

    for chunk_start in range(0, matrix_len, chunk_size):
        pairwise_similarity_chunk = similarity_by_chunk(chunk_start, chunk_start + chunk_size, X_test, X_train, matrix_len)
    
    dm.append(pairwise_similarity_chunk)
    Dm = pd.DataFrame(pairwise_similarity_chunk).round(2)
    sim = 1 - Dm.T
    sim = sim.reset_index()   
    
    descriptor_space = sim.rename(columns={'index': 'nearest_neighbor', 0:'similarity'})
    descriptor_space = descriptor_space.join(moldf[['InChIKey','Mol']])
    
    return descriptor_space

# Integrate Organ Weight for the top X nearest neighbors

def clustering(moldf, morgan_space, maccs_space, avalon_space, atom_pair_space):
    
    # Filter out compounds with no biological data
    morgan_drug_matrix = pd.merge(moldf[['InChIKey', 'DrugMatrix', 'compound_name']].dropna(axis=0, subset=['DrugMatrix']), morgan_space, how='inner', on='InChIKey')
    morgan_tg_gates = pd.merge(moldf[['InChIKey', 'TG-GATEs', 'compound_name']].dropna(axis=0, subset=['TG-GATEs']), morgan_space, how='inner', on='InChIKey')
    
    avalon_drug_matrix = pd.merge(moldf[['InChIKey', 'DrugMatrix', 'compound_name']].dropna(axis=0, subset=['DrugMatrix']), avalon_space, how='inner', on='InChIKey')
    avalon_tg_gates = pd.merge(moldf[['InChIKey', 'TG-GATEs', 'compound_name']].dropna(axis=0, subset=['TG-GATEs']), avalon_space, how='inner', on='InChIKey')
    
    # Filter out compounds with no biological data
    maccs_drug_matrix = pd.merge(moldf[['InChIKey', 'DrugMatrix', 'compound_name']].dropna(axis=0, subset=['DrugMatrix']), maccs_space, how='inner', on='InChIKey')
    maccs_tg_gates = pd.merge(moldf[['InChIKey', 'TG-GATEs', 'compound_name']].dropna(axis=0, subset=['TG-GATEs']), maccs_space, how='inner', on='InChIKey')
    
    # Filter out compounds with no biological data
    atom_pair_drug_matrix = pd.merge(moldf[['InChIKey', 'DrugMatrix', 'compound_name']].dropna(axis=0, subset=['DrugMatrix']), atom_pair_space, how='inner', on='InChIKey')
    atom_pair_tg_gates = pd.merge(moldf[['InChIKey', 'TG-GATEs', 'compound_name']].dropna(axis=0, subset=['TG-GATEs']), atom_pair_space, how='inner', on='InChIKey')
       
    # Merge clustered chemicals with biological data from each source
    morgan_drug_matrix_nns = pd.merge(morgan_drug_matrix, drug_matrix, on='DrugMatrix')
    morgan_tg_gates_nns = pd.merge(morgan_tg_gates, tg_gates, on='TG-GATEs')

    maccs_drug_matrix_nns = pd.merge(maccs_drug_matrix, drug_matrix, on='DrugMatrix')
    maccs_tg_gates_nns = pd.merge(maccs_tg_gates, tg_gates, on='TG-GATEs')
    
    atom_pair_drug_matrix_nns = pd.merge(atom_pair_drug_matrix, drug_matrix, on='DrugMatrix')
    atom_pair_tg_gates_nns = pd.merge(atom_pair_tg_gates, tg_gates, on='TG-GATEs')
        
    avalon_drug_matrix_nns = pd.merge(avalon_drug_matrix, drug_matrix, on='DrugMatrix')
    avalon_tg_gates_nns = pd.merge(avalon_tg_gates, tg_gates, on='TG-GATEs')
    
    
    # Preparing data for visualization
    
    morgan_nns_df = pd.concat([morgan_drug_matrix_nns, morgan_tg_gates_nns], axis=0).sort_values(by='similarity', ascending=False)
    morgan_nns_df = morgan_nns_df.join(pd.DataFrame({'descriptor': ['morgan'] * morgan_nns_df.shape[0]}))
      
    maccs_nns_df = pd.concat([maccs_drug_matrix_nns, maccs_tg_gates_nns], axis=0).sort_values(by='similarity', ascending=False)
    maccs_nns_df = maccs_nns_df.join(pd.DataFrame({'descriptor': ['maccs'] * maccs_nns_df.shape[0]}))
    
    atom_pair_nns_df = pd.concat([atom_pair_drug_matrix_nns, atom_pair_tg_gates_nns], axis=0).sort_values(by='similarity', ascending=False)
    atom_pair_nns_df = atom_pair_nns_df.join(pd.DataFrame({'descriptor': ['atom_pair'] * atom_pair_nns_df.shape[0]}))
    
    avalon_nns_df = pd.concat([avalon_drug_matrix_nns, avalon_tg_gates_nns], axis=0).sort_values(by='similarity', ascending=False)
    avalon_nns_df = avalon_nns_df.join(pd.DataFrame({'descriptor': ['avalon'] * avalon_nns_df.shape[0]}))
    
    all_desc = pd.concat([morgan_nns_df, maccs_nns_df, avalon_nns_df, atom_pair_nns_df], axis=0)
    available_descriptors = all_desc['descriptor'].unique()
    
    df = all_desc.groupby(['InChIKey',
                           'compound_name',
                           'parameter_type',
                           'parameter',
                           'dose',
                           'time',
                           'time_unit',
                           'value_unit',
                           'descriptor',
                           'source'], as_index=False).agg({'value': 'median',
                                                           'outcome': 'median',
                                                           'similarity': 'median'
                                                          })
    
    df = df.round(2)
       
    return df

# Check if metal
def metal_atomic_numbers(at):
    n = at.GetAtomicNum()
    return (n==13) or (n>=21 and n<=31) or (n>=39 and n<=50) or (n>=57 and n<=83) or (n>=89 and n<=115)

def is_metal(smiles):
    mol = Chem.MolFromSmiles(smiles)
    rwmol = Chem.RWMol(mol)
    rwmol.UpdatePropertyCache(strict=False)
    metal = [at for at in rwmol.GetAtoms() if metal_atomic_numbers(at)]
    return metal

# Get mol from SMILES input
def get_mol(smiles):
    if smiles is None:
        raise Exception('Invalid or empty input.')
        
    smiles = smiles.replace('@','')

    mol0 = Chem.MolFromSmiles(smiles)
    
    if mol0 is None:
        raise Exception('Invalid or empty input.')
    
    if is_metal(smiles):
        raise Exception('The inserted structure contains metal(s).')
           
    mol1 = standardizer.standardize_mol(mol0)
    smol, _ = standardizer.get_parent_mol(mol1)
    
    if '.' in Chem.MolToSmiles(smol):
        raise Exception('The input may be a mixture. Please insert one molecule at a time.')

    return smol

def get_similarity_df(mol):   
    # Calculate descriptors for molecule
    morgan, maccs, atom_pair, avalon  = calculate_descriptors(mol)

    # Get the nearest neighbors in all the descriptor spaces
    morgan_space = get_nearest_neighbor(morgan_train, morgan)
    maccs_space = get_nearest_neighbor(maccs_train, maccs)
    atom_pair_space = get_nearest_neighbor(atom_pair_train, atom_pair)
    avalon_space = get_nearest_neighbor(avalon_train, avalon)

    # Get the clustering dataset
    df = clustering(moldf, morgan_space, maccs_space, atom_pair_space, avalon_space)
    
    return df