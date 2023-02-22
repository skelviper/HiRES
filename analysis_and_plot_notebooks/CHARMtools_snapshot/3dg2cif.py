# transform 3dg/xyz to mmcif
# @Date 210723

import sys
import pandas as pd

def cli(args):
    input_file, output_file, factorBpath, maxGap= \
        args.input_file, args.output_file, args.factorBpath, args.maxGap
    threedg_to_cif(input_file, output_file, factorBpath, maxGap)

class Bin:
    def __init__(self,id:int,chromosome:str,position:int,x_coord:float,y_coord:float,z_coord:float,factorB:float=None):
        self.id = id
        self.chromosome = chromosome
        self.position = position
        self.nextBin = None
        self.previousBin = None
        self.x_coord = x_coord
        self.y_coord = y_coord
        self.z_coord = z_coord
        self.factorB = factorB
   
    def outputAtomLine(self,atomNum:int):
        """
        """
        return "\t".join(["HETATM",".",str(atomNum),self.chromosome,"003",
                          "1",str(self.position),str(self.x_coord),str(self.y_coord),str(self.z_coord),
                          str(self.factorB if self.factorB else "."),self.chromosome])
    def outputBondLine(self,bondNum:int,maxGap:int):
        """
        return bond of this bin and it's next when they are connected
        """
        if (self.chromosome == self.nextBin.chromosome and abs(self.nextBin.position-self.position) <= maxGap):
            return "\t".join([str(bondNum),"covale",self.chromosome,"003",
                          "1",str(self.position),self.nextBin.chromosome,"003",
                          "1",str(self.nextBin.position)])
        else: return None
        
def threedg_to_cif(tdgPath:str,outputCifPath:str,factorBpath:str=None,maxGap:int=1000000):
    """
    funtion as its name
    """
    positions = pd.read_table(tdgPath,header=None,names="chr pos x y z".split()).replace("[()]","",regex=True) # read in
    if(factorBpath!=None):
        factorB = pd.read_table(factorBpath,header=None,names="chrom pos factorB".split())
        positions = pd.merge(positions.assign(chrom=positions.chr.replace('.at','',regex=True)),factorB,how="left")
        print(positions)
    grouped = positions.groupby("chr") # split chromosomes
    
    binList = []
    binNum = 0
    for chr_name, coord_per_chrom in grouped:
        for index, series in coord_per_chrom.iterrows():
            # create atom for each bin
            if (factorBpath):
                currentBin = Bin(index,chr_name,series['pos'],series["x"],series["y"],series["z"],series["factorB"])
            else : currentBin = Bin(index,chr_name,series['pos'],series["x"],series["y"],series["z"])
            if(binNum !=0):
                binList[-1].nextBin = currentBin
            binNum += 1
            # if b factor is specificed.
            binList.append(currentBin)
    #print(binList)
    binList[-1].nextBin = binList[0]
    
    with open(outputCifPath,"w") as output_cif:
        output_cif.write("data_"+output_cif.name.replace(".cif","")+"\n")
        output_cif.write("#\n")
        output_cif.write("#\n")
        #write all connnection
        output_cif.write("""loop_\n_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
""")
        bondIndex = 1
        for bin in binList:
            #print(bondIndex)
            temp = bin.outputBondLine(bondIndex,maxGap)
            if temp:
                output_cif.write(temp + "\n")
                bondIndex += 1

        #write all atoms
        output_cif.write("""##
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.id
_atom_site.label_asym_id
_atom_site.label_comp_id
_atom_site.label_seq_id
_atom_site.label_atom_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.B_iso_or_equiv
_atom_site.auth_asym_id
""")
        atomIndex = 1
        for bin in binList:
            output_cif.write(bin.outputAtomLine(atomIndex)+"\n")
            atomIndex += 1


if __name__ == '__main__':
    threedg_to_cif(sys.argv[1],sys.argv[2],sys.argv[3])
    #for debugging
    #threedg_to_cif("/share/home/zliu/shareb/zliu/analysis/hires_gastrulation/HiC/post_m_structures/example/GasbE751168.10k.ma.chr1m.3dg",
    #    "/share/home/zliu/shareb/zliu/analysis/hires_gastrulation/HiC/post_m_structures/example/GasbE751168.10k.ma.chr1m.cif")
