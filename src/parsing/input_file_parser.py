import types

'''
[ForceFeild]
protein_ff = leaprc.protein.ff14SB
small_molecule_ff = leaprc.gaff
lipid_ff = leaprc.lipid17
water_ff = leaprc.water.tip3p
externel_amber_parms = "zn-mg.frcmod TPO.frcmod SO4.frcmod"
externel_amber_prep = "zn.prepi oh.prepi TPO.prepi SO4.prepi"
[Externel_program_path]
amber = $AMBERHOME
gaussian = $g03root
amber_md = pmemd.cuda_SPFP
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
ifcomplex_prepare = True
charge_model = AM1-BCC #or RESP, which need gaussian
tleap_file = None
ifoptimization_ligand = False
iffull_auto = True
'''

class InputParser:
    def __init__(self, filename):
        # Define the expected keys for each section and their default values
        self.def_sections()

        # Initialize variables to hold the parsed data
        self.data = {}
        for section, options in self.sections.items():
            self.data[section] = {}
            for option, default_value in options.items():
                self.data[section][option] = default_value

        self.current_section = None

        # Parse the input file
        with open(filename, 'r') as f:
            for line in f:
                # Remove inline comments
                line = line.split('#')[0].strip()
                if not line:  # Skip empty lines and comments
                    continue
                if line.startswith('[') and line.endswith(']'):
                    self.current_section = line[1:-1]
                elif self.current_section is not None:
                    for item in line.split(','):
                        key, value = item.split('=')
                        key = key.strip()
                        value = value.strip()
                        if key in self.sections[self.current_section]:
                            if value.isdigit():
                                self.data[self.current_section][key] = int(value)
                            elif value.replace('.', '').isdigit():
                                self.data[self.current_section][key] = float(value)
                            elif value.lower() == 'true':
                                self.data[self.current_section][key] = True
                            elif value.lower() == 'false':
                                self.data[self.current_section][key] = False
                            elif value.lower() == 'none':
                                self.data[self.current_section][key] = None
                            else:
                                self.data[self.current_section][key] = value.strip("'").strip('"')


        # Define the section-related methods dynamically
        for section in self.sections:
            def get_section_data(self, section=section):
                return self.data.get(section, {})
            setattr(self, f'get_{section}', types.MethodType(get_section_data, self))

    def def_sections(self):
        # Define the expected keys for each section and their default values
        self.sections = {
                            'ForceField': {
                                'protein_ff': 'leaprc.protein.ff14SB',
                                'small_molecule_ff': 'leaprc.gaff',
                                'lipid_ff': 'leaprc.lipid17',
                                'water_ff': 'leaprc.water.tip3p',
                                'externel_amber_parms': 'zn-mg.frcmod TPO.frcmod SO4.frcmod',
                                'externel_amber_prep': 'zn.prepi oh.prepi TPO.prepi SO4.prepi'
                            },
                            'Externel_program_path': {
                                'automd_home': '/home/zli/automd',
                                'amber': '$AMBERHOME',
                                'gaussian': '$g03root',
                                'amber_md': 'pmemd.cuda_SPFP',
                                'gaussian_scr_path': '/tmp/zli/scr'
                            },
                            'MD_input_parameters': {
                                'restraint_lig': '',
                                'restraint_rec': '',
                                'restraint_add': '',
                                'cutoff': 10.0,
                                'final_temp': 300.0,
                                'heat_nstlim': 25000,
                                'density_nstlim': 25000,
                                'equil_nstlim': 50000,
                                'prod_nstlim': 1000000,
                                'prod_steps': 4,
                            },
                            'Preparation_option': {
                                'ligand_pdb': 'mol.pdb',
                                'lig_residue_name': 'MOL',
                                'ligand_charge': 0,
                                'receptor_pdb': 'protein.pdb',
                                'ifonly_small_molecule_prepare': False,
                                'ifonly_protein_prepare': False,
                                'ifcomplex_prepare': True,
                                'charge_model': 'bcc',
                                'boxtype': 'solvateoct',
                                'boxsize': 10.0,
                                'sbond_file': None,
                                'tleap_file': None,
                                'ifoptimization_ligand': False,
                                'iffull_auto': True,
                            },
                            'Non_standard_amino_acid_prep_option': {
                                'ifonly_nonstandard_aminoacid_prepare': False,
                                'nonstandard_aminoacid_pdb': 'nonstandard_aminoacid.pdb',
                                'nonstandard_aminoacid_resname': 'UNK',
                                'nonstandard_aminoacid_residx': 1,
                                'nonstandard_aminoacid_charge_model': 'resp',
                                'nonstandard_aminoacid_charge': 0,
                                'ifoptimization_nonstandard_aminoacid': False,
                            },
        }



if __name__ == '__main__':
    def main():
        filename = 'input.txt'
        parser = InputParser(filename)
        ForceFeild_option = parser.get_ForceField()
        print(ForceFeild_option)
        Externel_program_path = parser.get_Externel_program_path()
        print(Externel_program_path)
        MD_input_parameters = parser.get_MD_input_parameters()
        print(MD_input_parameters) 
        Preparation_option = parser.get_Preparation_option()    
        print(Preparation_option)
    main()

