import sys
import os
if len(sys.argv) > 1:
    pdbfile=sys.argv[1]
    resname=sys.argv[2]



class Atom():
    def __init__(self,atmname,resname,x,y,z,chg=0.0,seq=0,resid=1):
        self.atmname=atmname
        self.resname=resname
        self.coord_x=float(x)
        self.coord_y=float(y)
        self.coord_z=float(z)
        self.chg=float(chg)
        self.seq=int(seq)
        self.resid=int(resid)
        self.special_list=['CL','Cl','BR','Br']
        if atmname[0:2] in self.special_list:
            self.atmtype=atmname[0:2]
        else:
            self.atmtype=atmname[0]
        
    def set_xyz(self, coord_x, coord_y, coord_z):
        self.coord_x = float(coord_x)
        self.coord_y = float(coord_y)
        self.coord_z = float(coord_z)

    def set_seq(self, seq):
        self.seq = int(seq)

    def set_chg(self, chg):
        self.chg = float(chg)
        
    def calc_dist(self, atom2):
        dist_sqr = (self.coord_x - atom2.coord_x)**2 + (self.coord_y - atom2.coord_y)**2 + (self.coord_z - atom2.coord_z)**2
        return dist_sqr ** 0.5
        
    def __repr__(self):
        return "Atom('" + ':%-2s@%-2s' % (self.resname,self.atmname) + "')"
    def __str__(self):
        return "ATOM  %5d  %-3s%4s%6d%12.3f%8.3f%8.3f%10.6f" % (self.seq,self.atmname,self.resname,self.resid,self.coord_x,self.coord_y,self.coord_z,self.chg)

class AtomChg():
    def __init__(self,atmname,chg):
        self.atmname=atmname
        self.chg=float(chg)
        
def printambmask(atmlst):
    out=":"+atmlst[0].resname+"@"
    for atom in atmlst:
        out+=atom.atmname
        out+=','
    print(out)

def format_lig_pdb(pdbfile, resname):
    if len(resname) == 1:
        resname='LG'+resname
    elif len(resname) == 2:
        resname='L'+resname
    elif len(resname) >= 4:
        resname=resname[0:3]

    f=open(pdbfile)
    addid=1
    namelst=[]
    lineend_lst=[]
    special_list=['CL','Cl','BR','Br']
    for line in f:
        line_end=line.strip()[26:]
        line_lst=line.strip().split()
        if (line_lst[0] == 'ATOM' or line_lst[0] == 'HETATM') and (len(line_lst) >= 8):
            #deal with duplicate names
            if addid >= 10: 
                print("Too many duplicated names")
                sys.exit(1)
            name=line_lst[2]
            if name[0:2] in special_list:
                atmtype=name[0:2]
                backup_name=atmtype+str(addid)
            else:
                atmtype=name[0]
                backup_name=atmtype+'X'+str(addid)
            if name in namelst:
                name=backup_name
                addid+=1

            namelst+=[name]
            lineend_lst+=[line_end]

    fo=open('ligand-format.pdb','w')
    i=1
    for name in namelst:
        print("ATOM  %5d  %-3s %3s%6d%s" % (i,name,resname,1,lineend_lst[i-1]),file=fo)
        i+=1

def pdb_to_xyz(ligfile, xyzfile):
    f=open(ligfile)
    LG_lst=[]
    for line in f:
        line=line.strip().split()
        if (line[0] == 'ATOM' or line[0] == 'HETATM') and (len(line) >= 8) :
            LG_lst+=[Atom(line[2],line[3],line[5],line[6],line[7])]
    f.close()
    fo=open(xyzfile,'w')
    xyz_content_list = []
    for atom in LG_lst:
        print(" %-2s %10.5f%10.5f%10.5f" % (atom.atmtype,atom.coord_x,atom.coord_y,atom.coord_z),file=fo)
        xyz_content_list.append(" %-2s %10.5f%10.5f%10.5f\n" % (atom.atmtype,atom.coord_x,atom.coord_y,atom.coord_z))
    print("", file=fo)
    xyz_content_list.append("\n")
    return xyz_content_list

def gen_ligpdb_by_prepi_formated_pdb(formated_pdb, prepi_file, ligpdb):
    f=open(formated_pdb)
    LG_dir={}
    for line in f:
        line=line.strip().split()
        if (line[0] == 'ATOM' or line[0] == 'HETATM') and (len(line) >= 8) :
            LG_dir[line[2]]=Atom(line[2],line[3],line[5],line[6],line[7])
    f.close()

    fp=open(prepi_file)
    chg_lst=[]
    for line in fp:
        line=line.strip().split()
        if len(line) == 11 and line[1] != 'DUMM':
            chg_lst+=[AtomChg(line[1],line[-1])]
    fp.close()

    fo=open(ligpdb, 'w')
    i=0
    for atmchg in chg_lst:
        i=i+1
        LG_dir[atmchg.atmname].set_seq(i)
        LG_dir[atmchg.atmname].set_chg(atmchg.chg)
        print(LG_dir[atmchg.atmname], file=fo)

    