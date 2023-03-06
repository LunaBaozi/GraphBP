import os, copy, csv
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

# AURORA KINASE B KNOWN INHIBITORS following https://www.guidetopharmacology.org/GRAC/ObjectDisplayForward?objectId=1937
# interessante la presenza del fluoro in molti inibitori
# interessante anche come scompaia la @ nella forma canonica delle smiles


def addHs2molfromsmiles(smiles, filename='Hmol.sdf', output_type=None):  

    # Generate a 3D structure from smiles
    mol = Chem.MolFromSmiles(smiles)
    hmol = Chem.AddHs(mol)
    writer = Chem.SDWriter(filename)  
    writer.write(hmol)
    writer.close()

    if output_type == 'smiles':
        return Chem.MolToSmiles(hmol)
    elif output_type == 'mol':
        return hmol
    #elif output_type == 'sdf'
    return  


# ciau = pd.read_csv('/home/luna/Documents/Coding/GraphBP/GraphBP/aurkb_inhibitors.csv', header=0)
# ciau2 = pd.read_csv('/home/luna/Documents/Coding/GraphBP/GraphBP/aurkb_inhibitors.csv', header=0)
# ciao = dict(zip(known_aurkb_inhibitors.inhibitor, known_aurkb_inhibitors.smiles))
# aurkb_inhibitors = copy.deepcopy(ciao)
# aurkb_inhibitors = {key:[] for key in aurkb_inhibitors}


# print(known_aurkb_inhibitors.iloc[:, 0])    # target id
# print(known_aurkb_inhibitors.iloc[:, 1])    # target uniprot
# print(known_aurkb_inhibitors.iloc[:, 2])    # target species
# print(known_aurkb_inhibitors.iloc[:, 3])    # ligand
# print(known_aurkb_inhibitors.iloc[:, 4])    # ligand id
# print(known_aurkb_inhibitors.iloc[:, 5])    # ligand species
# print(known_aurkb_inhibitors.iloc[:, 6])    # ligand pubchem cid
# print(known_aurkb_inhibitors.iloc[:, 7])    # smiles

known_aurkb_inhibitors = pd.read_csv('/home/luna/Documents/Coding/GraphBP/GraphBP/aurdata/aurorakinaseBinteractions.csv', header=0)
inhibitors_dict = dict(zip(known_aurkb_inhibitors.iloc[:, 3], known_aurkb_inhibitors.iloc[:, 7]))
# correspondence_dict = copy.deepcopy(inhibitors_dict)
# correspondence_dict = {key:[] for key in correspondence_dict}


def tanimoto_calc(smi1, smi2):
    
    mol1 = Chem.MolFromSmiles(smi1)
    mol2 = Chem.MolFromSmiles(smi2)
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048)

    s = round(DataStructs.TanimotoSimilarity(fp1, fp2), 3)

    return s



def tanimoto_calculator_from_smiles():

    tanimoto_canon = [['knowninhib_name', 'knowninhib_smiles', 'genmol_name', 'genmol_smiles', 'simil']]

    gen_mols_dir = '../trained_model/aurkb_gen_complex_1k'

    for gen_mol_dir in os.listdir(gen_mols_dir):
        rel_path = os.path.join(gen_mols_dir, gen_mol_dir)
        
        for filename in os.listdir(rel_path):
            if 'ligand' in filename:
                sppl = Chem.SDMolSupplier(os.path.join(rel_path, filename)) 
                
                for mol in sppl:
                    if mol is not None:  # some compounds cannot be loaded.
                        basic_smiles_pattern = Chem.MolToSmiles(mol)
                        canon_smiles_pattern = Chem.CanonSmiles(basic_smiles_pattern)
                        
                        for name, compound in inhibitors_dict.items():
                            # print(compound)
                            basic_smiles_test = str(compound)
                            canon_smiles_test = Chem.CanonSmiles(compound)

                            # using basic smiles
                            simil = tanimoto_calc(basic_smiles_pattern, basic_smiles_test)
                            tanimoto_basic.append([name, basic_smiles_test, filename, basic_smiles_pattern, simil])

                            # using canon smiles
                            simil = tanimoto_calc(canon_smiles_pattern, canon_smiles_test)
                            tanimoto_canon.append([name, canon_smiles_test, filename, canon_smiles_pattern, simil])

    return tanimoto_basic, tanimoto_canon 




tanimoto_basic, tanimoto_canon = compare_mol_smiles()

with open('./aurdata/results/output_tanimoto_basic_1k.csv', 'w', newline='') as csvfile1:
    writer = csv.writer(csvfile1)
    writer.writerows(tanimoto_basic)

with open('./aurdata/results/output_tanimoto_canon_1k.csv', 'w', newline='') as csvfile2:
    writer = csv.writer(csvfile2)
    writer.writerows(tanimoto_canon)

