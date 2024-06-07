import os
import sys
from ..nonstandard_residue_preparation import antechamber_relate_module
from ..nonstandard_residue_preparation import gaussian_relate_module
from ..nonstandard_residue_preparation import formate_lig_pdb


def main(forcefield_needed, input_lig_pdb, lig_resname, charge_model='bcc', lig_net_charge=0, gaussian_scr_path=None, gaussian_excute=None, ligifopt=False,
            ):
    """
    The workflow of the solvated receptor-ligand system preparation for amber MD simulation.

    :param forcefield_needed: dictioanry, the forcefield needed for the system.
    :param input_lig_pdb: string, the input ligand pdb file.
    :param lig_resname: string, the ligand residue name.
    :param charge_model: string, the charge model, 'bcc' or 'resp'.
    :param lig_net_charge: int, the net charge of the ligand.
    :param gaussian_scr_path: string, the path of the gaussian scratch.
    :param gaussian_excute: string, the excute file of the gaussian.
    :param ligifopt: bool, whether to optimize the ligand structure.
    """
    protein_ff = forcefield_needed['protein_ff']
    # Generate the ligand parameters
    small_molecule_ff = forcefield_needed['small_molecule_ff']
    formate_lig_pdb.format_lig_pdb(input_lig_pdb, lig_resname)# Generate the ligand-format.pdb
    ifcheck_gaussian_excute = True if charge_model == 'resp' else False
    gen_gaussian = gaussian_relate_module.Format_pdb_gen_gaussian('ligand-format.pdb', 1, lig_resname, lig_net_charge, protein_ff, gaussian_scr_path, gaussian_excute, ifcheck_gaussian_excute)# Generate the MOL_qm_gaussian.pdb
    gaussian_lig_pdb = f"{lig_resname}_qm_gaussian.pdb"
    if charge_model == 'resp':
        gen_gaussian.create_gaussian_com(ifopt=ligifopt)# Generate the gaussian input file
        gaussian_input_name = gen_gaussian.gaussian_input_name
        chk_file = gen_gaussian.chk_file
        run_gaussian = gaussian_relate_module.Gaussian_run(gaussian_input_name, gaussian_excute, chk_file)
        run_gaussian.run_gaussian()# Run the gaussian
        gaussian_out = f"{lig_resname}_qm.log"
    elif charge_model == 'bcc':
        gaussian_lig_pdb = f"{lig_resname}_qm_gaussian.pdb"
        gaussian_out = f"{lig_resname}_qm.log"    
    prepi_generator = antechamber_relate_module.PrepiGenerator(gaussian_lig_pdb, lig_resname, lig_net_charge, charge_model, small_molecule_ff.split('.')[-1], gaussian_out)
    prepi_generator.gen_small_molecule_prepi() 
    if os.path.exists(f"{lig_resname}.prepi"):
        if os.path.exists(f"{lig_resname}.frcmod"):
            print("Ligand parameters successfully generated.")
    else:
        print("Ligand parameters generation failed.")
        sys.exit()
  