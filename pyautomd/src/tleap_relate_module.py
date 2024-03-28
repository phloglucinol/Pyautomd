import os
import subprocess
import shutil
from .executable_checker import Executable_checker

# TODO: add the copy operator to the class to copy the externel parameters files.

class Tleap_runner:
    def __init__(self, automd_home, forcefields, external_amberparms, external_amberprep, rec_name_u, lig_name, boxtype, com_boxsize, lig_boxsize, sbond_file):
        self.automd_home = automd_home
        self.forcefields = forcefields # "leaprc.protein.ff14SB leaprc.water.tip3p leaprc.gaff leaprc.lipid17".split()
        self.external_amberparms = external_amberparms # "zn-mg.frcmod TPO.frcmod SO4.frcmod".split()
        self.external_amberprep = external_amberprep # "zn-mg.prepi TPO.prepi SO4.prepi".split()
        self.rec_name_u = rec_name_u # "protein"
        self.lig_name = lig_name # "MOL"
        self.boxtype = boxtype # "solvateoct" [solvatebox|solvateoct]
        self.com_boxsize = com_boxsize
        self.lig_boxsize = lig_boxsize
        self.sbond_file = sbond_file # "sbond.lst_bypdb_preparer"
        self.use_programs = ["tleap"]
        self.check_programs()
        self.copy_external_parameters()

    def copy_external_parameters(self):
        """Copies external parameter and prep files to the current directory."""
        source_directory = os.path.join(self.automd_home, "src", "external_parameters")

        # Combine parameter and prep files into a single list
        all_external_files = self.external_amberparms + self.external_amberprep

        # Copy files from the combined list
        for file in all_external_files:
            source_path = os.path.join(source_directory, file)
            if os.path.isfile(source_path):
                shutil.copy(source_path, '.')
            else:
                print(f"Warning: File {file} not found in {source_directory}")

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

    def run_tleap(self, auto):
        if not auto:
            self.run_tleap_direct()
        else:
            self.generate_tleap_file()
            self.run_tleap_direct()

    def run_tleap_direct(self):
        print("Running tleap, please wait...")
        # Check if tleap.txt exists
        if os.path.isfile("tleap.txt"):
            # Run tleap if the file exists
            subprocess.run(["tleap", "-f", "tleap.txt"], stdout=open("tleap.log", "w"))
            # Check if output file exists and is not empty
            if os.path.isfile(f"{self.rec_name_u}.prmtop"):
                print("tleap success")
            else:
                print("tleap error")
                exit(0)
        else:
            # Print error message if tleap.txt does not exist
            print("Error: tleap.txt does not exist.")
            exit(0)


    def generate_tleap_file(self):
        print("Running protein tleap, please wait...")
        with open("tleap.txt", "w") as file:
            for parm in self.forcefields:
                file.write(f"source {parm}\n")
            file.write(f"loadamberparams {self.lig_name}.frcmod\n")
            file.write(f"loadamberprep {self.lig_name}.prepi\n")
            for parm in self.external_amberparms:
                file.write(f"loadamberparams {parm}\n")
            for prep in self.external_amberprep:
                file.write(f"loadamberprep {prep}\n")

        # Replace placeholders in template and append to tleap.txt
        with open(os.path.join(self.automd_home, 'src', 'template_files', 'tleap.template'), "r") as template_file, open("tleap.txt", "a") as file:
            for line in template_file:
                line = line.replace("_rec_name_u_", self.rec_name_u).replace("_rec_name_l_", self.rec_name_u)
                line = line.replace("_BOXTYPE_", self.boxtype).replace("_COM_BOXSIZE_", str(self.com_boxsize))
                line = line.replace("_BOXTYPE_", self.boxtype).replace("_LIG_BOXSIZE_", str(self.lig_boxsize))
                file.write(line)


        if self.sbond_file is not None:
            with open(self.sbond_file, "r") as file:
                sbond_line = file.read()
        # Placeholder replacement for sbonds
        with open("tleap.txt", "r") as file:
            lines = file.readlines()
        with open("tleap.txt", "w") as file:
            for line in lines:
                if "_SBOND_LINE_" in line:
                    if self.sbond_file is not None:
                        file.write(sbond_line)
                    else:
                        file.write("\n")
                else:
                    file.write(line)
        

# Example usage
# tleap_runner = TLeapRunner(automd_home="/path/to/AUTOMD", forcefields=["forcefield1", "forcefield2"], external_amberparms=["parm1", "parm2"], external_amberprep=["prep1", "prep2"], rec_name_u="REC_U", rec_name_l="REC_L", lig_name="LIG", boxtype="BOX", boxsize="SIZE", sbond_file="sbond.txt")
# tleap_runner.run_tleap(auto="1")
