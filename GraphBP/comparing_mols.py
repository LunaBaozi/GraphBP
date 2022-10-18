import os, copy, csv
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

# AURORA KINASE B KNOWN INHIBITORS following https://www.guidetopharmacology.org/GRAC/ObjectDisplayForward?objectId=1937
# interessante la presenza del fluoro in molti inibitori
# interessante anche come scompaia la @ nella forma canonica delle smiles
# -> forse bisognerebbe provare a confrontare anche senza ridurre in forma canonica

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


# GREAT! THIS WORKS FOR COMPARING MOLECULES THAT HAVE DIFFERENT SMILES
# myPattern = 'CN1CCN(CC1)c1nc(Sc2ccc(cc2)NC(=O)C2CC2)nc(c1)Nc1[nH]nc(c1)C'
# myMolecule = '[H]c1c(N([H])c2c([H])c(C([H])([H])[H])nn2[H])nc(Sc2c([H])c([H])c(N([H])C(=O)C3([H])C([H])([H])C3([H])[H])c([H])c2[H])nc1N1C([H])([H])C([H])([H])N(C([H])([H])[H])C([H])([H])C1([H])[H]'
# a = Chem.CanonSmiles(myPattern)
# b = Chem.CanonSmiles(myMolecule)
# print(a)
# print(b)
# print(a==b)


# known_aurkb_inhibitors = pd.read_csv('/home/luna/Documents/Coding/GraphBP/GraphBP/aurkb_inhibitors.csv', header=0)
# # print(type(known_aurkb_inhibitors))  # pandas dataframe
# # known_aurkb_inhibitors.to_dict(orient='list')
# ciao = dict(zip(known_aurkb_inhibitors.inhibitor, known_aurkb_inhibitors.smiles))
# aurkb_inhibitors = copy.deepcopy(ciao)
# aurkb_inhibitors = {key:[] for key in aurkb_inhibitors}


known_aurkb_inhibitors = pd.read_csv('/home/luna/Documents/Coding/GraphBP/GraphBP/aurdata/aurorakinaseBinteractions.csv', header=0)
inhibitors_dict = dict(zip(known_aurkb_inhibitors.iloc[:, 3], known_aurkb_inhibitors.iloc[:, 7]))
correspondence_dict = copy.deepcopy(inhibitors_dict)
correspondence_dict = {key:[] for key in correspondence_dict}

# print(inhibitors_dict)
# print(list(known_aurkb_inhibitors.columns))
# print(known_aurkb_inhibitors.iloc[:, 0])    # target id
# print(known_aurkb_inhibitors.iloc[:, 1])    # target uniprot
# print(known_aurkb_inhibitors.iloc[:, 2])    # target species
# print(known_aurkb_inhibitors.iloc[:, 3])    # ligand
# print(known_aurkb_inhibitors.iloc[:, 4])    # ligand id
# print(known_aurkb_inhibitors.iloc[:, 5])    # ligand species
# print(known_aurkb_inhibitors.iloc[:, 6])    # ligand pubchem cid
# print(known_aurkb_inhibitors.iloc[:, 7])    # smiles
# print(known_aurkb_inhibitors.iloc[:, 8])
# print(known_aurkb_inhibitors.iloc[:, 9])
# print(known_aurkb_inhibitors.iloc[:, 10])
# print(known_aurkb_inhibitors.iloc[:, 11])



def tanimoto_calc(smi1, smi2):
    
    mol1 = Chem.MolFromSmiles(smi1)
    mol2 = Chem.MolFromSmiles(smi2)
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048)

    s = round(DataStructs.TanimotoSimilarity(fp1, fp2), 3)

    return s



def compare_mol_smiles():

    tanimoto = list()

    gen_mols_dir = '/home/luna/Documents/Coding/GraphBP/GraphBP/aurdata/aurkb_gen_complex'

    for gen_mol_dir in os.listdir(gen_mols_dir):
        rel_path = os.path.join(gen_mols_dir, gen_mol_dir)
        
        for filename in os.listdir(rel_path):
            if 'ligand' in filename:
                sppl = Chem.SDMolSupplier(os.path.join(rel_path, filename)) 
                
                for mol in sppl:
                    if mol is not None:  # some compounds cannot be loaded.
                        basic_smiles_pattern = Chem.MolToSmiles(mol)
                        # canon_smiles_pattern = Chem.CanonSmiles(basic_smiles_pattern)
                        
                        for name, compound in inhibitors_dict.items():
                            # print(compound)
                            basic_smiles_test = str(compound)
                            # canon_smiles_test = Chem.CanonSmiles(compound)

                            # should we use the standard or canonical smiles?
                            simil = tanimoto_calc(basic_smiles_pattern, basic_smiles_test)
                            tanimoto.append([name, basic_smiles_test, basic_smiles_pattern, simil])


                            # if basic_smiles_pattern == basic_smiles_test:
                            #     matches[compound].append(filename)
                            # elif basic_smiles_pattern == canon_smiles_test:
                            #     matches[compound].append(filename)
                            # elif canon_smiles_pattern == basic_smiles_test:
                            #     matches[compound].append(filename)
                            # elif canon_smiles_pattern == canon_smiles_test:
                            #     matches[compound].append(filename)
                            
                    # print('None')

    return tanimoto #matches




tanimoto = compare_mol_smiles()

# # for tlist in tanimoto:
# #     print(*tlist)

with open('output_withcsvfromwebsite.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(tanimoto)