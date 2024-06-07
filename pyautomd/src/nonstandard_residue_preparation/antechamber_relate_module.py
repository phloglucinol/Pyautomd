import os
import subprocess
import sys
import glob
import re
import shutil
from ..executable_checker import Executable_checker

class PrepiGenerator:
    """
    A class for generating prepi and frcmod files for nonstandard amino acids and small molecules.
    """
    def __init__(self, pdb, resn, netcharge, charge_method='bcc', lig_ff='gaff', gau=None):
        """
        Initializes the PrepiGenerator class.

        :param pdb: PDB file name. In pyautomd, usually be MOL_qm_gaussian.pdb 
        :param resn: Residue name. In pyautomd, usually be MOL
        :param charge_method: Charge method. 'bcc' or 'resp'
        :param lig_ff: Ligand force field. 'gaff' or 'gaff2'
        :param gau: Gaussian output file name. Required for resp charge method.
        :param netcharge: Net charge of the ligand. 
        """
        self.pdb = pdb
        self.resn = resn
        self.netcharge = str(int(netcharge))
        self.charge_method = charge_method # 'bcc' or 'resp'
        self.lig_ff = lig_ff # 'gaff' or 'gaff2'
        self.gau = gau
        if self.charge_method == 'resp':
            if self.gau is None:
                raise ValueError("Gaussian output file is required for resp charge method.")
        self.prefix = f"{resn}tmpwgXeY"
        # Define charge dictionaries
        self.ace = {"h1": 0.11230, "ch": -0.36620, "c": 0.59720, "o": -0.56790}
        self.nme = {"n": -0.41570, "h": 0.27190, "ch": -0.14900, "h1": 0.09760}
        self.ac = {"c": 0.59730, "o": -0.56790, "n": -0.41570, "h": 0.27190}
        self.use_programs = ["antechamber", "parmchk2", "prepgen", "resp"]
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

    def clean_up(self):
        """
        Removes specific files and file patterns as part of the cleanup process.
        """
        # File patterns to be removed
        file_patterns = [
            f"{self.prefix}.*", "ANTECHAMBER*", "NEWPDB.PDB",
            "PREP.INF", "punch", "qout", "QOUT", "esout", "ATOMTYPE.INF", 'sqm.*'
        ]

        for pattern in file_patterns:
            # Use glob to match the pattern
            for file_path in glob.glob(pattern):
                try:
                    os.remove(file_path)  # Remove the file
                except OSError as e:
                    pass
                #     print(f"Error while removing {file_path}: {e}")

    def gen_atominfo_template(self, input_file):
        """
        Returns a function that generates a template file by extracting lines starting with 'ATOM'.
        The input file is specified by the input_file parameter.
        """
        def generate_template():
            output_file = f"{self.prefix}.{input_file.split('.')[-1]}tmp"
            with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
                for line in infile:
                    if line.startswith("ATOM"):
                        outfile.write(line[:20] + '\n')
            return output_file

        return generate_template

    def run_antechamber(self, fi, fo, i, o, c, nc=None, rn=None, rf=None, at=None, cf=None):
        """
        Executes the antechamber command with the given options.

        :param fi: Input file format.
        :param fo: Output file format.
        :param i: Input file name.
        :param o: Output file name.
        :param c: Charge method.
        :param nc: (Optional) Net charge.
        :param rn: (Optional) Residue name.
        :param rf: (Optional) Residue topology file name.
        :param at: (Optional) Atom type.
        :param cf: (Optional) Charge file name.
        """
        # Base command
        command = ["antechamber", "-fi", fi, "-fo", fo, "-i", i, "-o", o, "-c", c]

        # Optional parameters
        if nc is not None:
            command.extend(["-nc", nc])
        if rn is not None:
            command.extend(["-rn", rn])
        if rf is not None:
            command.extend(["-rf", rf])
        if at is not None:
            command.extend(["-at", at])
        if cf is not None:
            command.extend(["-cf", cf])
        try:
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running antechamber: {e}")  

    def run_resp(self, i, o, e, t, q=None):
        """
        Executes the RESP command with the given options.

        :param i: Input file name for RESP.
        :param o: Output file name for RESP.
        :param e: espot
        :param t: qout.
        :param q: qin.
        """
        # Base command
        command = ["resp", "-O", "-i", i, "-o", o, "-e", e, "-t", t]

        # Optional parameter for charge file
        if q is not None:
            command.extend(["-q", q])

        try:
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running RESP: {e}")

    def has_nme_or_ace_cap(self):
        """
        Checks if the pdb file contains 'NME' or 'ACE' and returns the count.
        """
        cap_count = 0
        with open(self.pdb, 'r') as pdb_file:
            for line in pdb_file:
                if 'NME' in line or 'ACE' in line:
                    cap_count += 1
        return cap_count

    def run_parmchk2(self, i, f, o):
        """
        Executes the parmchk2 command with the given options.

        :param i: Input file name.
        :param f: Input file format.
        :param o: Output file name.
        """
        command = ["parmchk2", "-i", i, "-f", f, "-o", o]

        # Execute the command
        try:
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running parmchk2: {e}")

    def run_prepgen(self, i, o, f, m, rn, rf):
        """
        Executes the prepgen command with the given options.
        
        :param i: Input file name.
        :param o: Output file name.
        :param f: Input file format.
        :param m: mainchain file name.
        :param rn: Residue name.
        :param rf: Residue file name.
        """
        command = ["prepgen", "-i", i, "-o", o, "-f", f, "-m", m, "-rn", rn, "-rf", rf]

        # Execute the command
        try:
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running prepgen: {e}")

    def keep_atominfo_prepi_to_pdb(self, ifnonstandard_aminoa_acid=False):
        """
        Change the atom names in the .prepi file to match the atom names in the .pdb file.

        :param ifnonstandard_aminoa_acid: (Optional) If the residue is a nonstandard amino acid.
        """
        if ifnonstandard_aminoa_acid:
            with open(f"{self.prefix}.prepi", "r") as prepi_file:
                prepi_content = prepi_file.read()
            m1_count = prepi_content.count("+M")
            m2_count = prepi_content.count("-M")
        # Read contents of tmp1 and tmp2 files
        convert_name = {}
        with open(f"{self.prefix}.pdbtmp", "r") as file1, open(f"{self.prefix}.actmp", "r") as file2:
            for line1, line2 in zip(file1, file2):
                # Concatenate lines from both files
                combined_line = line1.strip() + line2.strip()
                
                # Extract residue name, pdb atom name, and ac atom name
                resname = combined_line[17:20]
                if resname != "ACE" and resname != "NME":
                    pdb_atm = combined_line[12:16].strip()
                    ac_atm = combined_line[33:37].strip()
                    convert_name[ac_atm] = pdb_atm
                    # print(convert_name)
        # Process the .prepi file
        output_lines = []

        with open(f"{self.prefix}.prepi", "r") as prepi_file:
            i = 0
            for line in prepi_file:
                i += 1
                fields = line.split()

                if len(line) > 60 and len(fields) == 11 and fields[2] != "DU":
                    if fields[1] not in convert_name:
                        raise KeyError(f"Key {fields[1]} not found in convert_name for line: {line}")
                    else:
                        new_line = f"{line[:6]}{convert_name[fields[1]]:<4}{line[10:]}"
                elif line.startswith("IMPROPER"):
                    new_line = line
                    i = 1000
                elif i > 1000 and i < 2000:
                    if len(fields) == 0:
                        new_line = "\n"
                        i = 0
                        if ifnonstandard_aminoa_acid:
                            if m1_count < 1:
                                output_lines.append(" CA   +M   C    O\n")
                            if m2_count < 1:
                                output_lines.append(" -M   CA   N    H\n")
                            output_lines.append("\n")
                            i = 0
                    else:
                        imprp = []
                        for atom in fields[:4]:
                            if "M" not in atom:
                                # If "M" is not in atom name, convert the atom name using convert_name dictionary
                                if atom not in convert_name:
                                    raise KeyError(f"Key {atom} not found in convert_name for line: {line}")
                                else:
                                    imprp.append(convert_name[atom])
                            else:
                                # If "M" is in atom name, use the atom name as it is
                                imprp.append(atom)
                        # imprp = [convert_name.get(atom, atom) for atom in fields[:4]]
                        new_line = f"{''.join(f'{atom:>5}' for atom in imprp)}\n"
                elif line.startswith("LOOP"):
                    new_line = line
                    i = 2000
                elif i > 2000:
                    if len(fields) == 0:
                        new_line = "\n"
                        i = 0
                    else:
                        if fields[0] not in convert_name:
                            raise KeyError(f"Key {fields[0]} not found in convert_name for line: {line}")
                        else:
                            new_line = f"{convert_name[fields[0]]:>5}{convert_name[fields[1]]:>5}\n"
                else:
                    new_line = line

                output_lines.append(new_line)
        # Apply sed replacements using regular expressions
        sed_replacements = [
            (r" 2HB ", " HB2 "), (r" 3HB ", " HB3 "), (r" 1H1 ", " H11 "), 
            (r" 2H1 ", " H12 "), (r" 3H1 ", " H13 "), (r" 1H2 ", " H21 "),
            (r" 2H2 ", " H22 "), (r" 3H2 ", " H23 "), (r" 1H5 ", " H51 "),
            (r" 2H5 ", " H52 "), (r" 3H5 ", " H53 "), (r" 1H1' ", " H1'1 "),
            (r" 2H1' ", " H1'2 "), (r" 3H1' ", " H1'3 "), (r" 1H2' ", " H2'1 "),
            (r" 2H2' ", " H2'2 "), (r" 3H2' ", " H2'3 "), (r" 1H5' ", " H5'1 "),
            (r" 2H5' ", " H5'2 "), (r" 3H5' ", " H5'3 "),
        ]

        processed_output = "".join(output_lines)
        for pattern, replacement in sed_replacements:
            processed_output = re.sub(pattern, replacement, processed_output)

        # Write the processed output to the .prepi file
        with open(f"{self.resn}.prepi", "w") as output_file:
            output_file.write(processed_output)

    def gen_nonstandard_amino_acid_prepi(self):
        """
        Generates the prepi file and frcmod file for nonstandard amino acids.
        """
        if self.charge_method == 'resp':
            #  1). run antechamber to get template ac file and resp input file
            self.run_antechamber(fi="gout", fo="ac", i=self.gau, o=f"{self.prefix}.ac", c="resp", nc=self.netcharge)
            if os.path.exists("ANTECHAMBER.ESP"):
                shutil.copy("ANTECHAMBER.ESP", f"{self.prefix}.esp")
            else:
                raise FileNotFoundError("ANTECHAMBER.ESP not found. Please check if the Gaussian output file is correct.")
            # Process ANTECHAMBER_RESP1.IN and ANTECHAMBER_RESP2.IN
            for resp_file, output_suffix in [("ANTECHAMBER_RESP1.IN", "step1.respin"), ("ANTECHAMBER_RESP2.IN", "step2.respin")]:
                with open(resp_file, "r") as infile, open(f"{self.prefix}.{output_suffix}", "w") as outfile:
                    go = False
                    for line in infile:
                        if line.strip() == "&end":
                            go = True
                            outfile.write(line)
                        elif go and len(line.strip()) == 0:
                            break
                        else:
                            outfile.write(line)
            #  2). add constraints to resp input files
            #      - constrain charges of NME and ACE exactly the same as in 
            #        AMBER database
            #      - constrain C-O and N-H in backbone the same as in AMBER
            #        database
            ace = self.ace
            nme = self.nme
            ac = self.ac
            with open(self.pdb, 'r') as pdb_file, open(f"{self.prefix}.tmp", 'w') as tmp_file:
                i = 0
                for line in pdb_file:
                    if line.startswith("ATOM"):
                        i += 1
                        resname = line[17:20].strip()
                        atmname = line[12:16].strip()

                    if resname == "ACE":
                        if "HH3" in line:
                            tmp_file.write(f"{1:5d}{ace['h1']:10.5f}\n{1:5d}{i:5d}\n")
                        if atmname == "C":
                            tmp_file.write(f"{1:5d}{ace['c']:10.5f}\n{1:5d}{i:5d}\n")
                        if atmname == "O":
                            tmp_file.write(f"{1:5d}{ace['o']:10.5f}\n{1:5d}{i:5d}\n")
                        if atmname == "CH3":
                            tmp_file.write(f"{1:5d}{ace['ch']:10.5f}\n{1:5d}{i:5d}\n")
                    elif resname == "NME":
                        if "HH3" in line:
                            tmp_file.write(f"{1:5d}{nme['h1']:10.5f}\n{1:5d}{i:5d}\n")
                        if atmname == "N":
                            tmp_file.write(f"{1:5d}{nme['n']:10.5f}\n{1:5d}{i:5d}\n")
                        if atmname == "H":
                            tmp_file.write(f"{1:5d}{nme['h']:10.5f}\n{1:5d}{i:5d}\n")
                        if atmname == "CH3":
                            tmp_file.write(f"{1:5d}{nme['ch']:10.5f}\n{1:5d}{i:5d}\n")
                    else:
                        if atmname == "N":
                            tmp_file.write(f"{1:5d}{ac['n']:10.5f}\n{1:5d}{i:5d}\n")
                        if atmname == "H":
                            tmp_file.write(f"{1:5d}{ac['h']:10.5f}\n{1:5d}{i:5d}\n")
                        if atmname == "C":
                            tmp_file.write(f"{1:5d}{ac['c']:10.5f}\n{1:5d}{i:5d}\n")
                        if atmname == "O":
                            tmp_file.write(f"{1:5d}{ac['o']:10.5f}\n{1:5d}{i:5d}\n")

            # Append additional lines
            with open(f"{self.prefix}.tmp", "a") as tmp_file:
                tmp_file.write("\n\n")


            # Append content to step1 and step2 respin files
            for step_file in [f"{self.prefix}.step1.respin", f"{self.prefix}.step2.respin"]:
                with open(f"{self.prefix}.tmp", 'r') as tmp_file, open(step_file, 'a') as resp_file:
                    shutil.copyfileobj(tmp_file, resp_file)
            self.run_resp(i=f"{self.prefix}.step1.respin", o=f"{self.prefix}.step1.respout", e=f"{self.prefix}.esp", t=f"{self.prefix}.step1.crg")
            self.run_resp(i=f"{self.prefix}.step2.respin", o=f"{self.prefix}.step2.respout", e=f"{self.prefix}.esp", q=f"{self.prefix}.step1.crg", t=f"{self.prefix}.step2.crg")
            self.run_antechamber(fi="gout", fo="ac", i=self.gau, o=f"{self.prefix}.ac.2", c="rc", cf=f"{self.prefix}.step2.crg", at="amber")
            os.remove(f"{self.prefix}.ac")
            os.rename(f"{self.prefix}.ac.2", f"{self.prefix}.ac")
            #  3). make mainchain file
            ac_atmname = {}
            with open(f"{self.prefix}.ac", "r") as ac_file:
                i = 0
                for line in ac_file:
                    if line.startswith("ATOM"):
                        i += 1
                        ac_atmname[i] = line[12:16].strip()

            # Process $pdb file and generate $prefix.mainchain
            with open(self.pdb, "r") as pdb_file, open(f"{self.prefix}.mainchain", "w") as mainchain_file:
                i = 0
                for line in pdb_file:
                    if line.startswith("ATOM"):
                        i += 1
                        resname = line[17:20].strip()
                        atmname = line[12:16].strip()

                        if resname not in ["ACE", "NME"]:
                            if atmname == "N":
                                mainchain_file.write(f"HEAD_NAME {ac_atmname[i]}\n")
                            elif atmname == "C":
                                mainchain_file.write(f"TAIL_NAME {ac_atmname[i]}\n")
                            elif atmname == "CA":
                                mainchain_file.write(f"MAIN_CHAIN {ac_atmname[i]}\n")
                        else:
                            mainchain_file.write(f"OMIT_NAME {ac_atmname[i]}\n")

            # Append predefined lines
            with open(f"{self.prefix}.mainchain", "a") as mainchain_file:
                mainchain_file.write("PRE_HEAD_TYPE C\n")
                mainchain_file.write("POST_TAIL_TYPE N\n")
                mainchain_file.write("CHARGE 0.0\n")
            #  4). remove capped residues and generage force field files
            # os.remove(f"{self.prefix}.prepi")
            self.run_prepgen(i=f"{self.prefix}.ac", o=f"{self.prefix}.prepi", f="prepi", m=f"{self.prefix}.mainchain", rn=self.resn, rf="")
        elif self.charge_method == 'bcc':
            raise NotImplementedError("Charge method 'bcc' is not supported for the parameters preparation nonstandard amino acid yet.")
        self.gen_atominfo_template(f"{self.prefix}.ac")()
        self.gen_atominfo_template(self.pdb)()
        self.keep_atominfo_prepi_to_pdb(ifnonstandard_aminoa_acid=True)
        self.run_parmchk2(f"{self.resn}.prepi", "prepi", f"{self.resn}.frcmod")
        self.clean_up()

    def gen_small_molecule_prepi(self):
        """
        Generates the prepi file for small molecules.
        """
        if self.charge_method == 'resp':
            self.run_antechamber(fi="gout", fo="ac", i=self.gau, o=f"{self.prefix}.ac", c="resp", nc=self.netcharge,  rn=self.resn, rf='\"\"',at=self.lig_ff)
            self.run_antechamber(fi="gout", fo="prepi", i=self.gau, o=f"{self.prefix}.prepi", c="resp", nc=self.netcharge, rn=self.resn, rf='\"\"', at=self.lig_ff)
        elif self.charge_method == 'bcc':
            self.run_antechamber(fi="pdb", fo="ac", i=self.pdb, o=f"{self.prefix}.ac", c="bcc", nc=self.netcharge, rn=self.resn, at=self.lig_ff)
            self.run_antechamber(fi="pdb", fo="prepi", i=self.pdb, o=f"{self.prefix}.prepi", c="bcc", nc=self.netcharge, rn=self.resn, at=self.lig_ff)
        else:
            raise ValueError("Invalid charge method.")
        self.gen_atominfo_template(f"{self.prefix}.ac")()
        self.gen_atominfo_template(self.pdb)()
        self.keep_atominfo_prepi_to_pdb()
        self.run_parmchk2(f"{self.resn}.prepi", "prepi", f"{self.resn}.frcmod")
        # self.clean_up()



if __name__ == "__main__":
    # # Test for small molecule (RESP)
    # pdb = 'MOL_qm.pdb'
    # resn = 'MOL'
    # charge_method='resp'
    # lig_ff='gaff2'
    # gau = 'MOL_qm.log'

    # # Test for small molecule (BCC)
    # pdb = 'MOL_qm_gaussian.pdb'
    # resn = 'MOL'
    # charge_method='bcc'
    # lig_ff='gaff2'
    # gau = None

    # Test for nonstandard amino acid
    pdb = 'LDM_qm_gaussian.pdb'
    resn = 'LDM'
    charge_method='resp'
    lig_ff=None
    gau = 'LDM_qm.log'
    prepi_generator = PrepiGenerator(pdb, resn, charge_method, lig_ff, gau)
    if prepi_generator.has_nme_or_ace_cap() > 0:
        prepi_generator.gen_nonstandard_amino_acid_prepi()
    else:
        prepi_generator.gen_small_molecule_prepi()
