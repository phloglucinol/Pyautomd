from .nonstandard_residue_preparation.formate_lig_pdb import gen_ligpdb_by_prepi_formated_pdb

class PDB_simple_processor:
    def __init__(self, prepi_file, ligand_format_pdb, receptor_pdb):
        self.prepi_file = prepi_file
        self.ligand_format_pdb = ligand_format_pdb
        self.receptor_pdb = receptor_pdb
        self.rec_name_u = receptor_pdb.split(".")[0]
        with open(self.receptor_pdb, "r") as file:
            lines = file.readlines()
        self.watline = next((i for i, line in enumerate(lines, 1) if ("HOH" in line or "WAT" in line) and ("HETATM" in line or "ATOM" in line)), None)

    def generate_lig_pdb(self):
        print("Generating lig.pdb")
        gen_ligpdb_by_prepi_formated_pdb(formated_pdb=self.ligand_format_pdb, prepi_file=self.prepi_file, ligpdb='lig.pdb')
        with open("lig.pdb", "a") as file:
            file.write("TER\n")

    def generate_rec_pdb(self):
        print("Generating rec.pdb")
        with open(self.receptor_pdb, "r") as file:
            lines = file.readlines()
        if self.watline:
            with open("rec.pdb", "w") as file:
                file.writelines(lines[:self.watline - 1])
        else:
            with open("rec.pdb", "w") as file:
                file.writelines(lines)

    def generate_rec_lig_pdb(self):
        print("Generating rec-lig.pdb")
        with open("rec-lig.pdb", "w") as outfile:
            with open("rec.pdb", "r") as infile:
                outfile.write(infile.read())
            with open("lig.pdb", "r") as infile:
                outfile.write(infile.read())

    def generate_combined_pdb(self):
        combined_pdb_name = f"{self.rec_name_u}_MOL.pdb"
        print(f"Generating {combined_pdb_name}")
        with open(combined_pdb_name, "w") as outfile:
            with open("rec.pdb", "r") as infile:
                outfile.write(infile.read())
            with open("lig.pdb", "r") as infile:
                outfile.write(infile.read())
            with open(self.receptor_pdb, "r") as infile:
                if self.watline:
                    outfile.writelines(infile.readlines()[self.watline:])

        # Removing lines starting with 'CONECT' and empty lines
        with open(combined_pdb_name, "r") as file:
            lines = file.readlines()
        with open(combined_pdb_name, "w") as file:
            file.writelines(line for line in lines if not line.startswith("CONECT") and line.strip())

