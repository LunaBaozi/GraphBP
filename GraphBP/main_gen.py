import pickle
from config import conf
from runner import Runner
import torch

runner = Runner(conf)


# POSSIAMO PROVARE A MODIFICARE QUESTOOOO
known_binding_site = True
# known_binding_site = False


node_temp = 0.5
dist_temp = 0.3
angle_temp = 0.4
torsion_temp = 1.0

min_atoms = 10
max_atoms = 45
focus_th = 0.5
contact_th = 0.5

# OCCHIO QUA SOTTO
num_gen = 20 #200 # number generate for each reference rec-lig pair

trained_model_path = 'trained_model'
epochs = [33]



for epoch in epochs:
    print('Epoch:', epoch)
    runner.model.load_state_dict(torch.load('{}/model_{}.pth'.format(trained_model_path, epoch)))
    all_mol_dicts = runner.generate(num_gen, 
                                    temperature=[node_temp, dist_temp, angle_temp, torsion_temp], 
                                    min_atoms=min_atoms,
                                    max_atoms=max_atoms,  
                                    focus_th=focus_th, 
                                    contact_th=contact_th, 
                                    add_final=True, 
                                    known_binding_site=known_binding_site)
    
    with open('{}/{}_mols.mol_dict'.format(trained_model_path, epoch),'wb') as f:
        pickle.dump(all_mol_dicts, f)
        


# NOTE: the structure of the data files is as follows:
# {class_label:d} {affinity:f} {RMSD:f} {receptor_path} {ligand_path}
