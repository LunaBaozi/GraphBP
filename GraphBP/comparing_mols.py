import os, copy
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

# TOZASERTIB
# tozasertib_c_addedH = addHs2molfromsmiles(tozasertib_c, 'tozasertib_c_addedH.sdf', 'smiles')
# print(tozasertib_c_addedH)

# print(converter('tozasertib_c_addedH.sdf'))
# tozasertib_i_addedH = converter(addHs2molfromsmiles(tozasertib_i, 'tozasertib_i_addedH.sdf', 'smiles'))
# # HESPERADIN
# hesperadin_c_addedH = converter(addHs2molfromsmiles(hesperadin_c, 'hesperadin_c_addedH.sdf', 'smiles'))
# hesperadin_i_addedH = converter(addHs2molfromsmiles(hesperadin_i, 'hesperadin_i_addedH.sdf', 'smiles'))
# # ZM447439
# zm447439_c_addedH = converter(addHs2molfromsmiles(zm447439_c, 'zm447439_c_addedH.sdf', 'smiles'))
# zm447439_i_addedH = converter(addHs2molfromsmiles(zm447439_i, 'zm447439_i_addedH.sdf', 'smiles'))



# from rdkit import Chem
# GREAT! THIS WORKS FOR COMPARING MOLECULES THAT HAVE DIFFERENT SMILES
# myPattern = 'CN1CCN(CC1)c1nc(Sc2ccc(cc2)NC(=O)C2CC2)nc(c1)Nc1[nH]nc(c1)C'
# myMolecule = '[H]c1c(N([H])c2c([H])c(C([H])([H])[H])nn2[H])nc(Sc2c([H])c([H])c(N([H])C(=O)C3([H])C([H])([H])C3([H])[H])c([H])c2[H])nc1N1C([H])([H])C([H])([H])N(C([H])([H])[H])C([H])([H])C1([H])[H]'

# a = Chem.CanonSmiles(myPattern)
# b = Chem.CanonSmiles(myMolecule)
# print(a)
# print(b)
# print(a==b)


known_aurkb_inhibitors = pd.read_csv('/home/luna/Documents/Coding/GraphBP/GraphBP/aurkb_inhibitors.csv', header=0)
# print(type(known_aurkb_inhibitors))  # pandas dataframe
# known_aurkb_inhibitors.to_dict(orient='list')
ciao = dict(zip(known_aurkb_inhibitors.inhibitor, known_aurkb_inhibitors.smiles))
aurkb_inhibitors = copy.deepcopy(ciao)
aurkb_inhibitors = {key:[] for key in aurkb_inhibitors}




def tanimoto_calc(mol1, mol2):   #(smi1, smi2):
    
    # mol1 = Chem.MolFromSmiles(smi1)
    # mol2 = Chem.MolFromSmiles(smi2)
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048)

    s = round(DataStructs.TanimotoSimilarity(fp1, fp2), 3)

    return s


# THERE MIGHT BE SOME NONE AROUND

def compare_mol_smiles():

    tanimoto = list()

    matches = copy.deepcopy(aurkb_inhibitors)
    # matches = {'tozasertib_c': [],
    #            'hesperadin_c': [],
    #            'zm447439_c': []
    #            }

    gen_mols_dir = '/home/luna/Documents/Coding/GraphBP/GraphBP/aurdata/aurkb_gen_complex'

    for gen_mol_dir in os.listdir(gen_mols_dir):
        rel_path = os.path.join(gen_mols_dir, gen_mol_dir)
        for filename in os.listdir(rel_path):
            if 'ligand' in filename:
                sppl = Chem.SDMolSupplier(os.path.join(rel_path, filename)) #if os.path.isfile(filename):
                for mol in sppl:
                    if mol is not None:  # some compounds cannot be loaded.
                        basic_smiles_pattern = Chem.MolToSmiles(mol)
                        canon_smiles_pattern = Chem.CanonSmiles(basic_smiles_pattern)
                        for compound in ciao.values():
                            # print(compound)
                            basic_smiles_test = compound
                            canon_smiles_test = Chem.CanonSmiles(compound)
                            # print(filename, pattern, test)

                            simil = tanimoto_calc(basic_smiles_pattern, basic_smiles_test)
                            tanimoto.append([(basic_smiles_test, basic_smiles_pattern, simil)])


                            if basic_smiles_pattern == basic_smiles_test:
                                matches[compound].append(filename)
                            elif basic_smiles_pattern == canon_smiles_test:
                                matches[compound].append(filename)
                            elif canon_smiles_pattern == basic_smiles_test:
                                matches[compound].append(filename)
                            elif canon_smiles_pattern == canon_smiles_test:
                                matches[compound].append(filename)
                            

    return matches


print(compare_mol_smiles())
# compare_mol_smiles()

