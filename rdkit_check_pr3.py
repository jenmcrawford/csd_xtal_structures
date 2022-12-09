import numpy as np
from scipy.spatial import distance
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import fnmatch
import os
import subprocess

# Chemistry
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from rdkit.Chem.Draw.MolDrawing import MolDrawing,DrawingOptions
DrawingOptions.bondLineWidth=1.8
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem import DataStructs

def pr3_match(pr3_ref,pr3_test):
    """
    Funtion to evaluate whether two molecules are identical through Tanimoto similarity
    arg1: smiles for reference molecule from UTS
    arg2: molecule to be evaluated from CSD

    returns: molecule to keep if each molecule is a substructure of the other
    """
    pr3_ref_mol = Chem.MolFromSmiles(pr3_ref)
    pr3_ref_Hs = Chem.AddHs(pr3_ref_mol)
    reffp = Chem.RDKFingerprint(pr3_ref_Hs)
    
    try:
        fps = Chem.RDKFingerprint(pr3_test)
        sim = DataStructs.FingerprintSimilarity(fps,reffp)
        if sim != 1.0:
            print(pr3_test.GetProp("_Name"))
            print(sim)
        elif sim == 1.0:
            return pr3_test
    except:
        with open('manual_check.txt','w') as f:
            f.write(pr3_test.GetProp("_Name"))

with open('phosphines_tosearch.txt') as f:
    full_list = f.readlines()
    gsk_id = [i.split(',',1)[0] for i in full_list]
    smi = [i.split(',',1)[1].strip() for i in full_list]

for i in range(len(gsk_id)):
    mol_list = []
    for file in os.listdir('.'):
	    if fnmatch.fnmatch(file,gsk_id[i] + '*.mol2'):
		    mol_list.append(file)
    pr3_csd = [Chem.MolFromMol2File(m, removeHs=False) for m in mol_list if m is not None]
    for j in range(len(pr3_csd)):
        if pr3_csd[j] is not None:
            pr3_csd[j].SetProp("_Name",mol_list[j].split('.')[0])
    matches = [pr3_match(smi[i],pr3_csd[j]) for j in range(len(pr3_csd)) if pr3_csd[j] is not None]
    filtered_matches = list(filter(None,matches))

    for m in filtered_matches:
        Chem.MolToMolFile(m, m.GetProp("_Name")+".mol")

     
    #Make image file of all matching phosphine ligands just for consistency check

    #Move files matching names with consistent SMILES to new folder

