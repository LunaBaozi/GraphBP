from rdkit import Chem


tozasertib_c = 'CN1CCN(CC1)c1nc(Sc2ccc(cc2)NC(=O)C2CC2)nc(c1)Nc1[nH]nc(c1)C'
tozasertib_i = 'CN1CCN(CC1)c1nc(Sc2ccc(cc2)NC(=O)C2CC2)nc(c1)Nc1[nH]nc(c1)C'
hesperadin_c = 'CCS(=O)(=O)Nc1ccc2c(c1)C(=C(c1ccccc1)Nc1ccc(cc1)CN1CCCCC1)C(=O)N2'
hesperadin_i = 'CCS(=O)(=O)Nc1ccc2c(c1)/C(=C(\c1ccccc1)/Nc1ccc(cc1)CN1CCCCC1)/C(=O)N2'
zm447439_c = 'COc1cc2c(ncnc2cc1OCCCN1CCOCC1)Nc1ccc(cc1)NC(=O)c1ccccc1'
zm447439_i = 'COc1cc2c(ncnc2cc1OCCCN1CCOCC1)Nc1ccc(cc1)NC(=O)c1ccccc1'


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
tozasertib_c_addedH = addHs2molfromsmiles(tozasertib_c, 'tozasertib_c_addedH.sdf', 'smiles')
print(tozasertib_c_addedH)

print(converter('tozasertib_c_addedH.sdf'))
# tozasertib_i_addedH = converter(addHs2molfromsmiles(tozasertib_i, 'tozasertib_i_addedH.sdf', 'smiles'))
# # HESPERADIN
# hesperadin_c_addedH = converter(addHs2molfromsmiles(hesperadin_c, 'hesperadin_c_addedH.sdf', 'smiles'))
# hesperadin_i_addedH = converter(addHs2molfromsmiles(hesperadin_i, 'hesperadin_i_addedH.sdf', 'smiles'))
# # ZM447439
# zm447439_c_addedH = converter(addHs2molfromsmiles(zm447439_c, 'zm447439_c_addedH.sdf', 'smiles'))
# zm447439_i_addedH = converter(addHs2molfromsmiles(zm447439_i, 'zm447439_i_addedH.sdf', 'smiles'))



# from rdkit import Chem
# GREAT! THIS WORKS FOR COMPARING MOLECULES THAT HAVE DIFFERENT SMILES
myPattern = 'CN1CCN(CC1)c1nc(Sc2ccc(cc2)NC(=O)C2CC2)nc(c1)Nc1[nH]nc(c1)C'
myMolecule = '[H]c1c(N([H])c2c([H])c(C([H])([H])[H])nn2[H])nc(Sc2c([H])c([H])c(N([H])C(=O)C3([H])C([H])([H])C3([H])[H])c([H])c2[H])nc1N1C([H])([H])C([H])([H])N(C([H])([H])[H])C([H])([H])C1([H])[H]'

a = Chem.CanonSmiles(myPattern)
b = Chem.CanonSmiles(myMolecule)

print(a)
print(b)

print(a==b)
