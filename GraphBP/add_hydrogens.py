import os
from rdkit import Chem

'/home/luna/Documents/Coding/GraphBP/GraphBP/aurdata/aurkb_gen_complex/4af3A1/Hmol.sdf'
write_dir = '/home/luna/Documents/Coding/GraphBP/GraphBP/aurdata/aurkb_gen_complex'

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


# noH_mol = '/home/luna/Documents/Coding/GraphBP/GraphBP/aurdata/gen_mols_epoch_33/content/GraphBP/GraphBP/trained_model/gen_mols_epoch_33/1.sdf'
# # converter(noH_mol)
# addHs2molfromsmiles(converter(noH_mol))