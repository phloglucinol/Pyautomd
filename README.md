# PyAutoMD: Python-Based Automatic Preparation for Amber Molecular Dynamics Input Files

![PyAutoMD Logo](https://github.com/phloglucinol/Pyautomd/blob/main/LOGO.png)

## 1. Important Notes
### 1.1 Preprocessing the Receptor Input PDB
It's crucial to preprocess the receptor input PDB file meticulously, which can be done using Maestro and pdb_preparer.py. This includes adding hydrogens, orienting polar hydrogens, determining protonation states of histidine (HIS), correctly naming capping atoms, and properly inserting a TER line at the end of each protein segment to indicate breaks in covalent bonds.

### 1.2 Key Directories in PyAutoMD
- `$PYAUTOMD/src/external_parameters`: Contains user-prepared parameter files for non-standard residues, including frcmod and prepi files. These files are copied to the execution directory by the program, referenced through "external_amber_parms" and "external_amber_prep" keywords in the automd_input.txt file.
- `$PYAUTOMD/src/template_files`: Holds template files for various purposes: tleap templates, Gaussian charge calculation, Amber MD control files (*.in), and shell script templates for submitting MD tasks under Amber.

## 2. Usage
### 2.1 General Usage of PyAutoMD
Detailed descriptions of the automd option are as follows:
```
pyautomd -i automd_input.txt
```

### 2.2 Automated Preparation of Topology and Coordinate Files of a Solvated Complex
For preparing the ligand-receptor complex system and its Amber MD input files, the automd_input.txt file should include:

```
[ForceField]
protein_ff = leaprc.protein.ff14SB # Force field for protein
small_molecule_ff = leaprc.gaff2 # Force field for small molecule
lipid_ff = leaprc.lipid17 # Force field for lipid
water_ff = leaprc.water.tip3p # Force field for water
external_amber_parms = "zn-mg.frcmod TPO.frcmod SO4.frcmod" # frcmod files in $PYAUTOMD/src/external_parameters
external_amber_prep = "zn.prepi oh.prepi TPO.prepi SO4.prepi" # prepi files in $PYAUTOMD/src/external_parameters

[External_Program_Path] 
gaussian = g03 # Gaussian execution program, needed for the resp charge model. Can be 'None' for bcc charge model.
amber_md = pmemd.cuda_SPFP # Amber MD execution program
gaussian_scr_path = /tmp/zli/scr # Scratch path for Gaussian

[MD_Input_Parameters]
restraint_lig = '' # Ligand ambmask for restraint (can be empty)
restraint_rec = '' # Receptor ambmask for restraint (can be empty)
restraint_add = '' # Additional ambmask for restraint (can be empty)
cutoff = 10.0 # Cutoff distance for nonbonded interactions
final_temp = 300.0 # Final temperature of MD simulation
heat_nstlim = 25000 # MD steps for heating process
density_nstlim = 25000 # MD steps for density equilibration
equil_nstlim = 50000 # MD steps for equilibration
prod_nstlim = 1000000 # MD steps for production
prod_steps = 4 # MD stages for production

[Preparation_Option]
ligand_pdb = mol.pdb # Ligand PDB file name
lig_residue_name = MOL # Ligand residue name
ligand_charge = 0 # Ligand net charge
receptor_pdb = protein.pdb # Receptor PDB file name
ifonly_small_molecule_prepare = False # Set True to prepare only the small molecule
ifonly_protein_prepare = False # Set True to prepare only the protein
ifcomplex_prepare = True # Set True to prepare the ligand-receptor complex
charge_model = bcc # or resp, requiring Gaussian
boxtype = 'solvateoct' # Type of simulation box (solvatebox or solvateoct)
boxsize = 10.0 # Size of simulation box
sbond_file = None # File name for disulfide bonds
tleap_file = tleap.txt # tleap script file name (default: tleap.txt)
ifoptimization_ligand = False # Set True to optimize the ligand before gaussian calculation
iffull_auto = True # Set True for a fully automatic process, generating tleap.txt automatically
```

### 2.3 Usage for the automatic preparation of the parameters files (prepi and frcmod) of the ligand
To prepare the parameters files (prepi and frcmod) of the ligand, the automd_input.txt file should include:

```
[ForceField]
protein_ff = leaprc.protein.ff14SB # Force field for protein
small_molecule_ff = leaprc.gaff2 # Force field for small molecule
lipid_ff = leaprc.lipid17 # Force field for lipid
water_ff = leaprc.water.tip3p # Force field for water
external_amber_parms = "zn-mg.frcmod TPO.frcmod SO4.frcmod" # frcmod files in $PYAUTOMD/src/external_parameters
external_amber_prep = "zn.prepi oh.prepi TPO.prepi SO4.prepi" # prepi files in $PYAUTOMD/src/external_parameters

[External_Program_Path] 
gaussian = g03 # Gaussian execution program, needed for the resp charge model. Can be 'None' for bcc charge model.
amber_md = pmemd.cuda_SPFP # Amber MD execution program
gaussian_scr_path = /tmp/zli/scr # Scratch path for Gaussian

[MD_Input_Parameters]
restraint_lig = '' # Ligand ambmask for restraint (can be empty)
restraint_rec = '' # Receptor ambmask for restraint (can be empty)
restraint_add = '' # Additional ambmask for restraint (can be empty)
cutoff = 10.0 # Cutoff distance for nonbonded interactions
final_temp = 300.0 # Final temperature of MD simulation
heat_nstlim = 25000 # MD steps for heating process
density_nstlim = 25000 # MD steps for density equilibration
equil_nstlim = 50000 # MD steps for equilibration
prod_nstlim = 1000000 # MD steps for production
prod_steps = 4 # MD stages for production

[Preparation_Option]
ligand_pdb = mol.pdb # Ligand PDB file name
lig_residue_name = MOL # Ligand residue name
ligand_charge = 0 # Ligand net charge
receptor_pdb = protein.pdb # Receptor PDB file name
ifonly_small_molecule_prepare = True # Set True to prepare only the small molecule
ifonly_protein_prepare = False # Set True to prepare only the protein
ifcomplex_prepare = False # Set True to prepare the ligand-receptor complex
charge_model = bcc # or resp, requiring Gaussian
boxtype = 'solvateoct' # Type of simulation box (solvatebox or solvateoct)
boxsize = 10.0 # Size of simulation box
sbond_file = None # File name for disulfide bonds
tleap_file = tleap.txt # tleap script file name (default: tleap.txt)
ifoptimization_ligand = False # Set True to optimize the ligand before gaussian calculation
iffull_auto = True # Set True for a fully automatic process, generating tleap.txt automatically
```

### 2.4 Usage for the automatic preparation of the parameters files (prepi and frcmod) of the nonstandard amino acid
To prepare the parameters files (prepi and frcmod) of the nonstandard amino acid, the automd_input.txt file should include:

```
[ForceField]
protein_ff = leaprc.protein.ff14SB # Force field for protein
small_molecule_ff = leaprc.gaff2 # Force field for small molecule
lipid_ff = leaprc.lipid17 # Force field for lipid
water_ff = leaprc.water.tip3p # Force field for water
external_amber_parms = "zn-mg.frcmod TPO.frcmod SO4.frcmod" # frcmod files in $PYAUTOMD/src/external_parameters
external_amber_prep = "zn.prepi oh.prepi TPO.prepi SO4.prepi" # prepi files in $PYAUTOMD/src/external_parameters

[Externel_program_path]
gaussian = g03
amber_md = pmemd.cuda_SPFP
gaussian_scr_path = /tmp/zli/scr

[MD_input_parameters]
restraint_lig = '' # ligand ambmask used for restraint this line can be empty
restraint_rec = '' # receptor ambmask used for restraint this line can be empty
restraint_add = '' # additional ambmask used for restraint this line can be empty
cutoff = 10.0
final_temp = 300.0
heat_nstlim = 25000
density_nstlim = 25000
equil_nstlim = 50000
prod_nstlim = 1000000
prod_steps = 4

[Preparation_option]
ligand_pdb = mol.pdb
lig_residue_name = MOL
ligand_charge = 0
receptor_pdb = protein.pdb
ifonly_small_molecule_prepare = False
ifonly_protein_prepare = False
ifcomplex_prepare = False
charge_model = bcc #or RESP, which need gaussian
boxtype = 'solvateoct' # [solvatebox|solvateoct]
boxsize = 10.0
sbond_file = None
tleap_file = None
ifoptimization_ligand = False
iffull_auto = True

[Non_standard_amino_acid_prep_option]
ifonly_nonstandard_aminoacid_prepare = True # Set True to prepare only the nonstandard amino acid
nonstandard_aminoacid_pdb = 'LDM_ini.pdb' # Nonstandard amino acid PDB file name
nonstandard_aminoacid_resname = 'LDM' # Nonstandard amino acid residue name
nonstandard_aminoacid_residx = 148 # Nonstandard amino acid residue index
nonstandard_aminoacid_charge_model = resp # charge model for nonstandard amino acid, which need gaussian, only resp is supported so far
nonstandard_aminoacid_charge = 1 # net charge for nonstandard amino acid
ifoptimization_nonstandard_aminoacid = False # Set True to optimize the nonstandard amino acid before gaussian calculation
```
