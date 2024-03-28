import subprocess
import os
# from formate_lig_pdb import pdb_to_xyz # for directly excute this script
from .formate_lig_pdb import pdb_to_xyz # for import this script
from ..executable_checker import Executable_checker


class Format_pdb_gen_gaussian:
    """
    Class to encapsulate the functionality of the formatgauss shell script.
    """
    def __init__(self, pdb_file='ligand-format.pdb', resid=1, resn='MOL', lig_net_charge=0, force_field_for_protein_terminal='leaprc.protein.ff14SB', gaussian_scr_path='/tmp/zli/scr', gaussian_excute='g03', ifcheck_gaussian_excute=True):
        self.pdb_file = pdb_file
        self.resid = resid
        self.resn = resn
        self.lig_net_charge = lig_net_charge
        self.force_field_for_protein_terminal = force_field_for_protein_terminal
        self.gaussian_scr_path = gaussian_scr_path
        self.nproc = os.cpu_count()
        self.prepare_pdb()
        self.generate_leap_input()
        if os.path.exists(f"{self.resn}_qm.pdb"):
            self.pdb_file = f"{self.resn}_qm.pdb"
        else:
            print(f"*Error: Failed to generate {self.resn}_qm.pdb")
            raise SystemExit
        self.gaussian_pdb = self.process_pdb_for_gaussian()
        self.gaussian_input_name = f"{os.path.splitext(self.pdb_file)[0]}.com"
        self.chk_file = os.path.join(self.gaussian_scr_path, f'{os.getpid()}.chk')
        self.gaussian_excute = gaussian_excute
        if not ifcheck_gaussian_excute:
            self.use_programs = ["tleap"]
        else:
            self.use_programs = ["tleap", self.gaussian_excute]
        self.check_programs()
    
    def check_programs(self):
        """
        Checks if the required programs are available in the system PATH and executable.
        """
        prog_ts = Executable_checker()
        for program in self.use_programs:
            if not prog_ts.is_executable_available(program):
                raise FileNotFoundError(f"{program} not found. Please check if the program is installed and available in the system PATH.")
            else:
                print(f"{program} is available in the system PATH and executable.")
        

    def prepare_pdb(self):
        """
        Prepare the PDB file by selecting specific residues and generating a temporary file.
        """
        with open(f'{self.pdb_file}.tmp', 'w') as tmp_file, open(self.pdb_file) as pdb:
            for line in pdb:
                if line.startswith(("ATOM", "HETATM")):
                    id = int(line[22:28].strip())
                    if id == self.resid:
                        tmp_file.write(line)
                    elif id == self.resid - 1:
                        if line[12:16].strip() == "CA":
                            tmp_file.write(f"{line[:12]} CH3 ACE{line[20:]}")
                        elif line[12:16].strip() == "C":
                            tmp_file.write(f"{line[:12]} C   ACE{line[20:]}")
                        elif line[12:16].strip() == "O":
                            tmp_file.write(f"{line[:12]} O   ACE{line[20:]}")
                    elif id == self.resid + 1:
                        if line[12:16].strip() == "CA":
                            tmp_file.write(f"{line[:12]} CH3 NME{line[20:]}")
                        elif line[12:16].strip() == "N":
                            tmp_file.write(f"{line[:12]} N   NME{line[20:]}")
                        elif line[12:16].strip() == "H":
                            tmp_file.write(f"{line[:12]} H   NME{line[20:]}")

    def generate_leap_input(self):
        """
        Generate leap input file and run tleap to generate the final PDB file (${LIG_RESN_NAME}_qm.pdb).
        """
        leap_input = f"""source {self.force_field_for_protein_terminal}
x = loadpdb {self.pdb_file}.tmp
savepdb x {self.resn}_qm.pdb 
quit"""

        with open(f'{self.pdb_file}.leap', 'w') as file:
            file.write(leap_input)
        
        subprocess.run(["tleap", "-f", f"{self.pdb_file}.leap"], stdout=subprocess.PIPE)

        os.remove(f"{self.pdb_file}.tmp")
        os.remove(f"{self.pdb_file}.leap")

    def get_self_defined_basis_set(self):
        """
        Get the self-defined basis set for the ligand.
        """
        gaussian_basis_set = [
            "-C -H -O -N -P -S -F -Cl -Br 0\n",
            "6-31G*\n",
            "****\n",
            "-I     0\n",
            "S   5   1.00\n",
            " 444750.0                    0.00089\n",
            "  66127.00                   0.00694\n",
            "  14815.00                   0.03609\n",
            "   4144.900                  0.13568\n",
            "   1361.200                  0.33878\n",
            "S   2   1.00\n",
            "    508.4400                 0.43659\n",
            "    209.5900                 0.18375\n",
            "S   1   1.00\n",
            "     81.959                  1.00000\n",
            "S   1   1.00\n",
            "     36.805                  1.00000\n",
            "S   1   1.00\n",
            "     13.495                  1.00000\n",
            "S   1   1.00\n",
            "      6.8859                 1.00000\n",
            "S   1   1.00\n",
            "      2.5520                 1.00000\n",
            "S   1   1.00\n",
            "      1.2088                 1.00000\n",
            "S   1   1.00\n",
            "      0.2734                 1.00000\n",
            "S   1   1.00\n",
            "      0.1009                 1.00000\n",
            "P   4   1.00\n",
            "   2953.600                  0.01221\n",
            "    712.6100                 0.08587\n",
            "    236.7100                 0.29493\n",
            "     92.63100                0.47849\n",
            "P   1   1.00\n",
            "     39.73200                1.00000\n",
            "P   1   1.00\n",
            "     17.27300                1.000000\n",
            "P   1   1.00\n",
            "      7.957000               1.000000\n",
            "P   1   1.00\n",
            "      3.152900               1.000000\n",
            "P   1   1.00\n",
            "      1.332800               1.000000\n",
            "P   1   1.00\n",
            "      0.494700               1.000000\n",
            "P   1   1.00\n",
            "      0.216000               1.000000\n",
            "P   1   1.00\n",
            "      0.082930               1.000000\n",
            "D   3   1.00\n",
            "    261.9500                 0.03144\n",
            "     76.73400                0.19028\n",
            "     27.55100                0.47247\n",
            "D   1   1.00\n",
            "     10.60600                1.000000\n",
            "D   1   1.00\n",
            "      3.421700               1.000000\n",
            "D   1   1.00\n",
            "      1.137000               1.000000\n",
            "D   1   1.00\n",
            "      0.302000               1.000000\n",
            "****\n",
            "\n",
            "I 1.98\n", 
            "\n"
        ]
        return gaussian_basis_set

    def process_pdb_for_gaussian(self):
        """
        Processes the PDB file for Gaussian input.
        It extracts ATOM and HETATM lines, converts HETATM to ATOM,
        checks for NME/ACE caps, and identifies duplicate atom names.
        """
        gaussian_pdb = f'{os.path.splitext(self.pdb_file)[0]}_gaussian.pdb'

        # Extract ATOM and HETATM lines and convert HETATM to ATOM
        with open(self.pdb_file, 'r') as infile, open(gaussian_pdb, 'w') as outfile:
            for line in infile:
                if line.startswith(('ATOM', 'HETATM')):
                    outfile.write(line.replace('HETATM', 'ATOM  '))

        # Check for NME/ACE caps
        with open(gaussian_pdb, 'r') as file:
            cap_count = sum(1 for line in file if 'NME' in line or 'ACE' in line)

        if cap_count > 0:
            print("* Found NME/ACE Cap!")

        # Check for duplicate atom names excluding NME and ACE
        atom_names = []
        with open(gaussian_pdb, 'r') as file:
            for line in file:
                if 'NME' not in line and 'ACE' not in line:
                    atom_names.append(line[12:16].strip())

        duplicates = set([name for name in atom_names if atom_names.count(name) > 1])
        if duplicates:
            print(f"*Error: There are {len(duplicates)} duplicated atom names in \"{self.pdb_file}\"")
            print("*They are:", ' '.join(duplicates))
            raise SystemExit
        return gaussian_pdb  # Successful processing

    def convert_pdb_to_xyz(self, pdb_file):
        """
        Convert the PDB file to a Gaussian XYZ format.
        """
        xyz_content_list = pdb_to_xyz(pdb_file, f'{os.path.splitext(pdb_file)[0]}.xyz')
        return xyz_content_list

    def create_gaussian_com(self, ifopt=False):
        """
        Create Gaussian .com file based on the shell script details.
        """
        outname = self.gaussian_input_name
        charge = self.lig_net_charge
        nproc = self.nproc
        if ifopt == True:
            opt_str = 'opt'
        else:
            opt_str = 'sp'
        chk_file = self.chk_file
        com_content = [
            f'%nproc={nproc}\n',
            '%mem=2Gb\n',
            f'%chk={chk_file}\n',
            f'# PM3MM {opt_str}\n\n',
            'test\n\n',
            f'{charge} 1\n'
        ]

        # Append molecular coordinates
        com_content.extend(self.convert_pdb_to_xyz(self.gaussian_pdb))

        # Append additional Gaussian parameters
        com_content.extend([
            '--link1--\n',
            f'%nproc={nproc}\n',
            '%mem=2Gb\n',
            f'%chk={chk_file}\n',
            '# HF/gen SCF=tight Pop=(MK,ReadRadii) iop(6/33=2,6/41=10,6/42=17,6/50=1) Geom=AllCheck\n\n',
        ])
        com_content.extend(self.get_self_defined_basis_set())

        with open(f'{outname}', 'w',) as com_file:
            com_file.writelines(com_content)

class Gaussian_run():
    '''
    Class to run Gaussian calculations and handle the output.
    '''
    def __init__(self, input_file, gaussian_excute, chk_file):
        self.gaussian_excute = gaussian_excute
        self.input_file = input_file
        self.output_file = f'{os.path.splitext(self.input_file)[0]}.log'
        self.chk_file = chk_file

    def run_gaussian(self):
        # Run Gaussian and generate output file
        subprocess.run([self.gaussian_excute, self.input_file], stdout=subprocess.PIPE)
        if self.check_normal_termination():
            print('Gaussian calculation completed successfully!')
        else:
            print('Gaussian calculation failed!')
            raise SystemExit
        self.remove_chk_file()

    def remove_chk_file(self):
        """
        Remove the checkpoint file.
        """
        try:
            os.remove(self.chk_file)
        except FileNotFoundError:
            pass       

    def check_normal_termination(self):
        """
        Checks if the last line of a file starts with 'Normal termination'.

        :param file_path: Path to the file to be checked.
        :return: True if the last line starts with 'Normal termination', False otherwise.
        """
        file_path = self.output_file
        try:
            with open(file_path, 'r') as file:
                for last_line in file:
                    pass
                return last_line.startswith(" Normal termination")
        except FileNotFoundError:
            print(f"Error: File not found - {file_path}")
            return False
        except Exception as e:
            print(f"Error: An error occurred - {e}")
            return False


if __name__ == '__main__':
    pdb_file='ligand-format.pdb'
    resid=1
    resn='MOL'
    lig_net_charge=0
    force_field_for_protein_terminal='leaprc.protein.ff14SB'
    gaussian_scr_path='/tmp/zli/scr'
    gen_gaussian = Format_pdb_gen_gaussian(pdb_file, resid, resn, lig_net_charge, force_field_for_protein_terminal, gaussian_scr_path)
    gen_gaussian.create_gaussian_com(True)
    gaussian_excute = 'g03'
    input_file = gen_gaussian.gaussian_input_name
    chk_file = gen_gaussian.chk_file
    run_gaussian = Gaussian_run(input_file, gaussian_excute, chk_file)
    run_gaussian.run_gaussian()