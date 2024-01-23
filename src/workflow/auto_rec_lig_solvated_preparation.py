import os
import sys
from ..tleap_relate_module import Tleap_runner
from ..md_input_generator import MD_input_prep
from ..rec_lig_com_pdb_generator import PDB_simple_processor
from ..nonstandard_residue_preparation import antechamber_relate_module
from ..nonstandard_residue_preparation import gaussian_relate_module
from ..nonstandard_residue_preparation import formate_lig_pdb


def main(automd_home, forcefield_needed, amber_md, external_amberparms, external_amberprep, boxtype, boxsize,
            input_lig_pdb, lig_resname, charge_model='bcc', iffull_auto=True, prod_steps=4, cutoff = 10.0,
            final_temp=300.0, heat_nstlim=1000, density_nstlim=1000, equil_nstlim=10000, prod_nstlim=10000,
            restraint_lig=None, restraint_rec=None, restraint_add=None,
            receptor_pdb='protein.pdb', sbond_file=None, 
            lig_net_charge=0, gaussian_scr_path=None, gaussian_excute=None, ligifopt=False,
            ):
    """
    The workflow of the solvated receptor-ligand system preparation for amber MD simulation.

    :param automd_home: string, the path of the automd package.
    :param forcefield_needed: dictioanry, the forcefield needed for the system.
    :param amber_md: string, the path of the amber MD package. like: pmemd.cuda
    :param external_amberparms: string, the path of the external amber frcmod files.
    :param external_amberprep: string, the path of the external amber prepi files.
    :param boxtype: string, the type of the box, 'solvateoct' or 'solvatebox'.
    :param boxsize: float, the distance between the solute and the edge of the box.
    :param input_lig_pdb: string, the input ligand pdb file.
    :param lig_resname: string, the ligand residue name.
    :param charge_model: string, the charge model, 'bcc' or 'resp'.
    :param iffull_auto: bool, whether to run the whole process automatically.
    :param prod_steps: int, the MD production times. like 4 for 4 segments of MD production.
    :param cutoff: float, the cutoff distance of the nonbonded interaction.
    :param final_temp: float, the final temperature of the MD simulation.
    :param heat_nstlim: int, the MD steps of the heating process.
    :param density_nstlim: int, the MD steps of the density equilibration process.
    :param equil_nstlim: int, the MD steps of the equilibration process.
    :param prod_nstlim: int, the MD steps of the production process.
    :param restraint_lig: string, the ligand residue number. used for the restraint ligand.
    :param restraint_rec: string, the receptor residue number range. used for the restraint receptor.
    :param restraint_add: string, the additional residue number range. used for the restraint additional residue.
    :param receptor_pdb: string, the receptor pdb file.
    :param sbond_file: string, the disulfide bond file.
    :param lig_net_charge: int, the net charge of the ligand.
    :param gaussian_scr_path: string, the path of the gaussian scratch.
    :param gaussian_excute: string, the excute file of the gaussian.
    :param ligifopt: bool, whether to optimize the ligand structure.
    """
    protein_ff = forcefield_needed['protein_ff']
    rec_name_u = receptor_pdb.split('.')[0]
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
    prepi_generator = antechamber_relate_module.PrepiGenerator(gaussian_lig_pdb, lig_resname, charge_model, small_molecule_ff.split('.')[-1], gaussian_out)
    prepi_generator.gen_small_molecule_prepi() 
    if os.path.exists(f"{lig_resname}.prepi"):
        if os.path.exists(f"{lig_resname}.frcmod"):
            print("Ligand parameters successfully generated.")
    else:
        print("Ligand parameters generation failed.")
        sys.exit()
    # Generate the ligand, receptor, receptor-ligand complex pdb file
    pdb_processor = PDB_simple_processor(f"{lig_resname}.prepi", "ligand-format.pdb", receptor_pdb)
    pdb_processor.generate_lig_pdb() # generate lig.pdb
    pdb_processor.generate_rec_pdb() # generate rec.pdb
    pdb_processor.generate_rec_lig_pdb() # generate rec-lig.pdb
    pdb_processor.generate_combined_pdb() # generate protein_MOL.pdb
    # Generate the tleap input file and run tleap
    tleap_runner = Tleap_runner(automd_home, list(forcefield_needed.values()), external_amberparms, external_amberprep,  rec_name_u, lig_resname, boxtype, boxsize, sbond_file)
    tleap_runner.run_tleap(iffull_auto)
    # Generate the MD input file
    md_input_generator = MD_input_prep(automd_home, amber_md, rec_name_u, lig_resname, prod_steps, cutoff, final_temp, heat_nstlim, density_nstlim, equil_nstlim, prod_nstlim, restraint_lig, restraint_rec, restraint_add)
    md_input_generator.run()
