import os 
import shutil 
from rdkit import Chem
# from add_hydrogens import converter
# from add_hydrogens import addHs2molfromsmiles


aur_protein_pdb = '4af3A'
rec_name = '4af3_A_rec'
aurkb_file = '/home/luna/Documents/Coding/GraphBP/GraphBP/aurdata/AURKB_HUMAN_54_344_0/4af3_A_rec.pdb'
gen_ligands_path = '/home/luna/Documents/Coding/GraphBP/GraphBP/aurdata/gen_mols_epoch_33/content/GraphBP/GraphBP/trained_model/gen_mols_epoch_33'
aurkb_gen_complex = '/home/luna/Documents/Coding/GraphBP/GraphBP/aurdata/aurkb_gen_complex'
write_dir = '/home/luna/Documents/Coding/GraphBP/GraphBP/aurdata/aurkb_gen_complex/molsHs_added'


os.makedirs(os.path.dirname(aurkb_gen_complex), exist_ok=True)
os.makedirs(os.path.dirname(write_dir), exist_ok=True)


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


def addHs2molfromsmiles(smiles):  

    # Generate a 3D structure from smiles
    mol = Chem.MolFromSmiles(smiles)
    hmol = Chem.AddHs(mol)
    writer = Chem.SDWriter(os.path.join(f'{write_dir}', 'Hmol.sdf'))
    writer.write(hmol)
    writer.close()

    return hmol 



for i in range(1, 101):
    os.makedirs(os.path.join(aurkb_gen_complex, f'{aur_protein_pdb}{i}'))

    subfolder = os.path.join(aurkb_gen_complex, f'{aur_protein_pdb}{i}')
    noH_mol = os.path.join(gen_ligands_path, f'{i}.sdf')

    shutil.copy(aurkb_file, subfolder)
    shutil.copy(f'{gen_ligands_path}/{i}.sdf', subfolder)
    # shutil.copy(f'{write_dir}/Hmol.sdf', os.path.join(aurkb_gen_complex, f'{aur_protein_pdb}{i}'))

    # Returns
    mol_addedHs = addHs2molfromsmiles(converter(noH_mol))
    shutil.copy(f'{write_dir}/Hmol_{i}.sdf', subfolder)

    os.rename(f'{subfolder}/{rec_name}.pdb', f'{subfolder}/{aur_protein_pdb}{i}_protein.pdb')
    # os.rename(f'{subfolder}/{i}.sdf', f'{subfolder}/{aur_protein_pdb}{i}_ligand.sdf')
    os.rename(f'{subfolder}/Hmol_{i}.sdf', f'{subfolder}/{aur_protein_pdb}{i}_ligand.sdf')

        
        
        
        
        
# os.makedirs(os.path.dirname(f'{aur_protein}_ligand{i}'), exist_ok=True)
# shutil.copy(gen_ligands_path, dest_fpath)s