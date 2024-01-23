from optparse import OptionParser
from src.parsing.input_file_parser import InputParser
from src.workflow import auto_rec_lig_solvated_preparation
from src.workflow import auto_lig_parm_preparation
from src.workflow import auto_nonstandard_aminoacid_parm_preparation

class optParser():  
    def __init__(self, fakeArgs):
        parser = OptionParser()
        parser.add_option('-i', '--input', dest='input', help="The file name of input, which recording the automd input option. Default: 'input.txt'", default='input.txt')
        if fakeArgs:
            self.option, self.args = parser.parse_args(fakeArgs)
        else:
            self.option, self.args = parser.parse_args()

class Ligand_parameter_generation():
    def __init__(self, ligand_pdb, lig_residue_name, ligand_charge, ligand_charge_model):
        self.ligand_pdb = ligand_pdb
        self.lig_residue_name = lig_residue_name
        self.ligand_charge = ligand_charge
        self.ligand_charge_model = ligand_charge_model

opts = optParser('')
parser = InputParser(opts.option.input)
ForceFeild_option = parser.get_ForceField()
protein_ff = ForceFeild_option['protein_ff']
small_molecule_ff = ForceFeild_option['small_molecule_ff']
lipid_ff = ForceFeild_option['lipid_ff']
water_ff = ForceFeild_option['water_ff']
externel_amber_parms = ForceFeild_option['externel_amber_parms']
externel_amber_prep = ForceFeild_option['externel_amber_prep']
Externel_program_path = parser.get_Externel_program_path()
automd_home = Externel_program_path['automd_home']
amber = Externel_program_path['amber']
gaussian = Externel_program_path['gaussian']
gaussian_scr_path = Externel_program_path['gaussian_scr_path']
amber_md = Externel_program_path['amber_md']
MD_input_parameters = parser.get_MD_input_parameters()
restraint_lig = MD_input_parameters['restraint_lig']
restraint_rec = MD_input_parameters['restraint_rec']
restraint_add = MD_input_parameters['restraint_add']
cutoff = MD_input_parameters['cutoff']
final_temp = MD_input_parameters['final_temp']
heat_nstlim = MD_input_parameters['heat_nstlim']
density_nstlim = MD_input_parameters['density_nstlim']
equil_nstlim = MD_input_parameters['equil_nstlim']
prod_nstlim = MD_input_parameters['prod_nstlim']
prod_steps = MD_input_parameters['prod_steps']
Preparation_option = parser.get_Preparation_option()
ligand_pdb = Preparation_option['ligand_pdb']
lig_residue_name = Preparation_option['lig_residue_name']
ligand_charge = Preparation_option['ligand_charge']
receptor_pdb = Preparation_option['receptor_pdb']
ifonly_small_molecule_prepare = Preparation_option['ifonly_small_molecule_prepare']
ifonly_protein_prepare = Preparation_option['ifonly_protein_prepare']
ifcomplex_prepare = Preparation_option['ifcomplex_prepare']
charge_model = Preparation_option['charge_model']
boxtype = Preparation_option['boxtype']
boxsize = Preparation_option['boxsize']
sbond_file = Preparation_option['sbond_file']
tleap_file = Preparation_option['tleap_file']
ifoptimization_ligand = Preparation_option['ifoptimization_ligand']
iffull_auto = Preparation_option['iffull_auto']
non_standard_aminoacid_option = parser.get_Non_standard_amino_acid_prep_option()
ifonly_nonstandard_aminoacid_prepare = non_standard_aminoacid_option['ifonly_nonstandard_aminoacid_prepare']
nonstandard_aminoacid_pdb = non_standard_aminoacid_option['nonstandard_aminoacid_pdb']
nonstandard_aminoacid_resname = non_standard_aminoacid_option['nonstandard_aminoacid_resname']
nonstandard_aminoacid_residx = non_standard_aminoacid_option['nonstandard_aminoacid_residx']
nonstandard_aminoacid_charge_model = non_standard_aminoacid_option['nonstandard_aminoacid_charge_model']
nonstandard_aminoacid_charge = non_standard_aminoacid_option['nonstandard_aminoacid_charge']
ifoptimization_nonstandard_aminoacid = non_standard_aminoacid_option['ifoptimization_nonstandard_aminoacid']



forcefield_needed = {
    "protein_ff": protein_ff,
    "small_molecule_ff": small_molecule_ff,
    "lipid_ff": lipid_ff,
    "water_ff": water_ff,
}
externel_amber_parms = externel_amber_parms.split()
externel_amber_prep = externel_amber_prep.split()


if __name__ == '__main__':
    if ifcomplex_prepare:
        auto_rec_lig_solvated_preparation.main(automd_home=automd_home,
                                                forcefield_needed=forcefield_needed,
                                                amber_md=amber_md,
                                                external_amberparms=externel_amber_parms,
                                                external_amberprep=externel_amber_prep,
                                                boxtype=boxtype,
                                                boxsize=boxsize,
                                                input_lig_pdb=ligand_pdb,
                                                lig_resname=lig_residue_name,
                                                charge_model=charge_model,
                                                iffull_auto=iffull_auto,
                                                prod_steps=prod_steps,
                                                cutoff=cutoff,
                                                final_temp=final_temp,
                                                heat_nstlim=heat_nstlim,
                                                density_nstlim=density_nstlim,
                                                equil_nstlim=equil_nstlim,
                                                prod_nstlim=prod_nstlim,
                                                restraint_lig=restraint_lig,
                                                restraint_rec=restraint_rec,
                                                restraint_add=restraint_add,
                                                receptor_pdb=receptor_pdb,
                                                sbond_file=sbond_file,
                                                lig_net_charge=ligand_charge,
                                                gaussian_scr_path=gaussian_scr_path,
                                                gaussian_excute=gaussian,
                                                ligifopt=ifoptimization_ligand,
        )
    elif ifonly_small_molecule_prepare:
        auto_lig_parm_preparation.main(forcefield_needed=forcefield_needed,
                                        input_lig_pdb=ligand_pdb,
                                        lig_resname=lig_residue_name,
                                        charge_model=charge_model,
                                        lig_net_charge=ligand_charge,
                                        gaussian_scr_path=gaussian_scr_path,
                                        gaussian_excute=gaussian,
                                        ligifopt=ifoptimization_ligand,
        )
    elif ifonly_nonstandard_aminoacid_prepare:
        auto_nonstandard_aminoacid_parm_preparation.main(protein_ff=protein_ff,
                                                            input_aminoacids_pdb=nonstandard_aminoacid_pdb,
                                                            resname=nonstandard_aminoacid_resname,
                                                            residx=nonstandard_aminoacid_residx,
                                                            charge_model=nonstandard_aminoacid_charge_model,
                                                            net_charge=nonstandard_aminoacid_charge,
                                                            gaussian_scr_path=gaussian_scr_path,
                                                            gaussian_excute=gaussian,
                                                            ifopt=ifoptimization_nonstandard_aminoacid,
        )

