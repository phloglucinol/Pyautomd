import os
import sys
from ..nonstandard_residue_preparation import antechamber_relate_module
from ..nonstandard_residue_preparation import gaussian_relate_module


def main(protein_ff, input_aminoacids_pdb, resname, residx, charge_model='resp', net_charge=0, gaussian_scr_path='/tmp/zli/scr', gaussian_excute='g03', ifopt=False,
            ):
    """
    The workflow of the solvated receptor-ligand system preparation for amber MD simulation.

    :param protein_ff: string, the forcefield of the protein.
    :param input_aminoacids_pdb: string, the input aminoacids pdb file. May include the neighboring residues of the target nonstandard amino acid.(at least including the bone atoms)
    :param resname: string, the nonstandard amino acid residue name.
    :param residx: int, the nonstandard amino acid residue index in the input_aminoacids_pdb file.
    :param charge_model: string, the charge model, so far only 'resp' is supported.
    :param net_charge: int, the net charge of the nonstandard amino acid.
    :param gaussian_scr_path: string, the path of the gaussian scratch.
    :param gaussian_excute: string, the excute file of the gaussian.
    :param ifopt: bool, whether to optimize the nonstandard amino acid structure.
    """
    gen_gaussian = gaussian_relate_module.Format_pdb_gen_gaussian(input_aminoacids_pdb, residx, resname, net_charge, protein_ff, gaussian_scr_path)# Generate the MOL_qm_gaussian.pdb
    gaussian_lig_pdb = f"{resname}_qm_gaussian.pdb"
    if charge_model == 'resp':
        gen_gaussian.create_gaussian_com(ifopt=ifopt)# Generate the gaussian input file
        gaussian_input_name = gen_gaussian.gaussian_input_name
        chk_file = gen_gaussian.chk_file
        run_gaussian = gaussian_relate_module.Gaussian_run(gaussian_input_name, gaussian_excute, chk_file)
        run_gaussian.run_gaussian()# Run the gaussian
        gaussian_out = f"{resname}_qm.log"
    elif charge_model == 'bcc':
        raise ValueError("The charge model 'bcc' is not supported for nonstandard amino acid so far.")
    prepi_generator = antechamber_relate_module.PrepiGenerator(gaussian_lig_pdb, resname, net_charge, charge_model, None, gaussian_out)
    if prepi_generator.has_nme_or_ace_cap() > 0:
        prepi_generator.gen_nonstandard_amino_acid_prepi()
    else:
        print("The nonstandard amino acid input pdb is not including  NME or ACE.")
        sys.exit()
    if os.path.exists(f"{resname}.prepi"):
        if os.path.exists(f"{resname}.frcmod"):
            print("Non-standard amino acid parameters successfully generated.")
    else:
        print("Non-standard amino acid parameters generation failed.")
        sys.exit()
