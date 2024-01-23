import os
import shutil
import subprocess

class MD_input_prep:
    def __init__(self, automd_home, exec_prog, rec_name_u, lig_resn_name, prod_steps, cutoff, final_temp, heat_nstlim, density_nstlim, equil_nstlim, prod_nstlim, restraint_lig, restraint_rec, restraint_add):
        self.automd_home = automd_home
        self.exec_prog = exec_prog
        self.rec_name_u = rec_name_u
        self.lig_resn_name = lig_resn_name
        self.submit_pbs_head = os.path.join(f"{self.automd_home}", "src", "template_files", "submit.template")
        self.prod_steps = prod_steps
        self.restraint_add = restraint_add
        self.cutoff = cutoff
        self.final_temp = final_temp
        self.heat_nstlim = heat_nstlim
        self.density_nstlim = density_nstlim
        self.equil_nstlim = equil_nstlim
        self.prod_nstlim = prod_nstlim
        self.restraint_lig = restraint_lig # the residue index of the ligand in the tleap-protein-wat.pdb
        self.restraint_rec = restraint_rec # the residue index range of the receptor in the tleap-protein-wat.pdb
        self.restraint_add = restraint_add # the residue index range of the additional restraint in the tleap-protein-wat.pdb


    def run(self):
        self.set_restraint_parameters()
        self.create_run_directory()
        self.generate_submit_pbs()
        self.generate_input_files()
        print("automd preparation finished")

    def set_restraint_parameters(self):
        if not self.restraint_lig:
            with open(f"tleap-{self.rec_name_u}-wat.pdb", "r") as file:
                for line in file:
                    if self.lig_resn_name in line:
                        self.restraint_lig = line.split()[4]
                        break
        if not self.restraint_rec:
            self.restraint_rec = f"1-{int(self.restraint_lig) - 1}" # assume the ligand is the next residue after the receptor

    def create_run_directory(self):
        os.makedirs("run", exist_ok=True)
        os.chdir("run")

    def generate_submit_pbs(self):
        print("Generating submit.pbs")
        with open("submit.pbs", "w") as file:
            with open(self.submit_pbs_head, "r") as head_file:
                file.write(head_file.read())
            with open(os.path.join(f"{self.automd_home}", "src", "template_files", "submit2.template"), "r") as template_file:
                for line in template_file:
                    line = line.replace("_EXEC_PROG_", self.exec_prog.replace("/", "\/"))
                    line = line.replace("_rec_name_l_", self.rec_name_u)
                    file.write(line)
            for i in range(1, self.prod_steps + 1):
                last_step = i - 1
                incoord = "equil.rst" if last_step == 0 else f"prod{last_step}.rst"
                file.write(f"{self.exec_prog} -O -i ./ins/prod.in -o prod{i}.out -p {self.rec_name_u}.prmtop -c {incoord} -r prod{i}.rst -x {self.rec_name_u}-prod{i}.mdcrd\n")
                file.write(f"gzip -9 {self.rec_name_u}-prod{i}.mdcrd&\n")
            file.write("wait\n")
        os.chmod("submit.pbs", 0o755)
        print("submit.pbs successfully generated")

    def generate_input_files(self):
        print("Generating MD ins files")
        shutil.copytree(os.path.join(self.automd_home, "src", "template_files", "ins"), os.path.join(".", "ins"), dirs_exist_ok=True)

        # Replace placeholders in ins files
        for filename in os.listdir("ins"):
            with open(f"ins/{filename}", "r") as file:
                content = file.read()

            if self.restraint_add == "":
                content = content.replace(",_RESTRAINT_ADD_", "")
            
            content = content.replace("_RESTRAINT_LIG_", self.restraint_lig)
            content = content.replace("_RESTRAINT_REC_", self.restraint_rec)
            content = content.replace("_RESTRAINT_ADD_", self.restraint_add)
            content = content.replace("_CUTOFF_", str(self.cutoff))
            content = content.replace("_FINAL_TEMP_", str(self.final_temp))

            # Special replacements for specific files
            if filename == "heat.in":
                content = content.replace("_HEAT_NSTLIM_", str(self.heat_nstlim))
            elif filename == "density.in":
                content = content.replace("_DENSITY_NSTLIM_", str(self.density_nstlim))
            elif filename == "equil.in":
                content = content.replace("_EQUIL_NSTLIM_", str(self.equil_nstlim))
            elif filename == "prod.in":
                content = content.replace("_PROD_NSTLIM_", str(self.prod_nstlim))

            with open(f"ins/{filename}", "w") as file:
                file.write(content)

        print("MD ins files successfully generated")
