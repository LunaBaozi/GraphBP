import os 
import shutil 
from rdkit import Chem


# aur_protein_pdb = '4af3A'
# rec_name = '4af3_A_rec'
# aurkb_file = '/home/luna/Documents/Coding/GraphBP/GraphBP/aurdata/AURKB_HUMAN_54_344_0/4af3_A_rec.pdb'
# gen_ligands_path = '/home/luna/Documents/Coding/GraphBP/GraphBP/aurdata/gen_mols_epoch_33_1k/content/GraphBP/GraphBP/trained_model/gen_mols_epoch_33/'
# aurkb_gen_complex = '/home/luna/Documents/Coding/GraphBP/GraphBP/aurdata/aurkb_gen_complex_1k'
# write_dir = '/home/luna/Documents/Coding/GraphBP/GraphBP/aurdata/aurkb_gen_complex_1k/molsHs_added'

aurkA_pdb = '4ztq'
aurkB_pdb = '4af3'
aurkA_rec_name = '4ztq_A_rec'
aurkB_rec_name = '4af3_A_rec'
aurkA_file_path = '../datav11/crossdock2020/AURKA_HUMAN_122_397_0/4ztq_A_rec.pdb'
aurkB_file_path = '../datav11/crossdock2020/AURKB_HUMAN_54_344_0/4af3_A_rec.pdb'
gen_ligands_path = '../trained_model/gen_mols_epoch_33/'
aurkA_gen_complexes = './results/AURKA_GEN_COMPLEX'
aurkB_gen_complexes = './results/AURKB_GEN_COMPLEX'

# for file in os.listdir(gen_ligands_path):
#     print(type(file))
#     print(file.endswith('.sdf'))
#     break
# # [os.remove(file) for file in os.listdir('/data/trainGBP/GraphBP/GraphBP/trained_model/gen_mols_epoch_33/') if file.endswith('.png')]
# [os.remove(file) for file in os.listdir(gen_ligands_path) if file.endswith('.png')]
# print(len(os.listdir(gen_ligands_path)))


def converter(file_name):
    sppl = Chem.SDMolSupplier(file_name)
    outname = file_name.replace(".sdf", ".txt")
    out_file = open(outname, "w")
    for mol in sppl:
        if mol is not None:  # some compounds cannot be loaded.
            smi = Chem.MolToSmiles(mol)
            name = mol.GetProp("_Name")
            out_file.write(f"{smi}\t{name}\n")

            out_file.close()
            return smi
    return None
        # else: break
    # out_file.close()
    # return smi


def addHs2molfromsmiles(smiles):  

    # Generate a 3D structure from smiles
    mol = Chem.MolFromSmiles(smiles)
    hmol = Chem.AddHs(mol)
    writer = Chem.SDWriter('Hmol.sdf')  #(os.path.join(f'{write_dir}', 'Hmol.sdf'))
    writer.write(hmol)
    writer.close()

    return hmol 


def createEquiBindDataset(aurora_type):

    os.makedirs(os.path.dirname(f'aurk{aurora_type}_gen_complexes'), exist_ok=True)
    # os.makedirs(os.path.dirnmae(aurkB_gen_complexes), exist_ok=True)

    for i in range(1, len(os.listdir(gen_ligands_path))-2):   # (-2) because the folder also contains two dictionaries
        os.makedirs(os.path.join(aurkA_gen_complexes, f'{aurkA_pdb}{i}'))
        os.makedirs(os.path.join(aurkB_gen_complexes, f'{aurkB_pdb}{i}'))

        subfolderA = os.path.join(aurkA_gen_complexes, f'{aurkA_pdb}{i}')
        subfolderB = os.path.join(aurkB_gen_complexes, f'{aurkB_pdb}{i}')
        noH_mol = os.path.join(gen_ligands_path, f'{i}.sdf')

        shutil.copy(aurkb_file, subfolder)
        shutil.copy(f'{gen_ligands_path}/{i}.sdf', subfolder)

        # Returns
        sppl = Chem.SDMolSupplier(noH_mol)
        outname = noH_mol.replace(".sdf", ".txt")
        out_file = open(outname, "w")
        for mol in sppl:
            if mol is not None:  # some compounds cannot be loaded.
                smi = Chem.MolToSmiles(mol)
                name = mol.GetProp("_Name")
                out_file.write(f"{smi}\t{name}\n")
            else: break

        out_file.close()
        # return smi
        # smi = converter(noH_mol)
        if smi is not None:
            mol_addedHs = addHs2molfromsmiles(smi)
            shutil.move('Hmol.sdf', subfolder)

            os.rename(f'{subfolder}/{rec_name}.pdb', f'{subfolder}/{aur_protein_pdb}{i}_protein.pdb')
            os.rename(f'{subfolder}/Hmol.sdf', f'{subfolder}/{aur_protein_pdb}{i}_ligand.sdf')

        
# createEquiBindDataset()