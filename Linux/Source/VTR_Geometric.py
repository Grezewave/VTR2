#V1.2
#Contacts match by geometric propertie, VMD(Vector Medium Distance)
#Some entry changes
#Write a plot .pml file
#Match skill improved and refined
#Most detailed analisys avaiable(by residue and color scale)
#Most Graphic detail and execution stats
#More run options(Chain filter)

import sys
import Classify
import Contacts
import OSfunct
import Plot
import time
import VTR_Functions as vtr
import PymolGen as py
import os

def main():
#parameters -------------------------------------------------------------------------------------------------------------------------
    start = time.time()
    protein1 = sys.argv[1]
    protein2 = sys.argv[2]
    _type = ''

    if (len(sys.argv) >= 9):
        if 'd' == sys.argv[8]:
            _type += '-d'   

    outname = "../Logs/"+protein1[protein1.rfind("/")+1:-4]+"x"+protein2[protein2.rfind("/")+1:-4]+_type+"Log.txt"

#build aligned protein and contacts list --------------------------------------------------------------------------------------------
    path = OSfunct.TMAlign(protein1,protein2)
    rtt_name = "../Data/" + path + protein1[protein1.rfind("/"):-4] + "_rotate.pdb"
    rtt_protein = Classify.classify(rtt_name)
    stc_protein = Classify.classify(protein2)
    rtt_contacts = Contacts.contacts(rtt_protein,protein1[protein1.rfind("/"):-4] + "_rotate",sys.argv[4],sys.argv[5])
    stc_contacts = Contacts.contacts(stc_protein,protein2[protein2.rfind("/"):-4],sys.argv[6],sys.argv[7])

#match contacts ---------------------------------------------------------------------------------------------------------------------
    matches,rtt_dismatches,stc_dismatches = vtr.match_contacts(rtt_contacts,stc_contacts,int(sys.argv[3]))
    end = time.time()

#write output -----------------------------------------------------------------------------------------------------------------------
    vtr.write_dismatch(protein1,protein2,rtt_dismatches,stc_dismatches,_type)
    out = open(outname, 'w')
    out.write("Match execution time: " + str(end-start)+"\n")
    result = str(len(matches)) + " matches found" + "\n"
    out.write(result)
    if (0 == len(matches)):
        out.write("\n")
        out.write("\n")
        out.write(protein1[protein1.rfind("/")+1:]+"\n")
        if "/" != sys.argv[4]:
            out.write(sys.argv[4]+"\n")
        else:
            out.write("All\n")
        if "/" != sys.argv[5]:
            out.write(sys.argv[5]+"\n")
        else:
            out.write("All\n")
        out.write(protein2[protein2.rfind("/")+1:]+"\n")
        if "/" != sys.argv[6]:
            out.write(sys.argv[6]+"\n")
        else:
            out.write("All\n")
        if "/" != sys.argv[7]:
            out.write(sys.argv[7]+"\n")
        else:
            out.write("All\n")
        out.write(sys.argv[3]+"\n")

        if "d" in _type:
            out.write("Detailed-")
        else:
            out.write("Simple  -")

        if "e" in _type:
            out.write("Equivalence\n")
        else:
            out.write("Distance\n")
        out.close()
        sys.exit()
    out.write("RMSD = "+str(vtr.RMSD(matches, rtt_protein, stc_protein)) + "\n")
    out.write("VTR = "+str(vtr.VTR(matches, rtt_protein, stc_protein,len(rtt_dismatches),len(stc_dismatches))) + "\n")
    vtr.writer(protein1,protein2,rtt_protein,stc_protein,rtt_contacts,stc_contacts,matches,_type)

#make plots -------------------------------------------------------------------------------------------------------------------------
    folder = py.detailed_ploter(rtt_name, protein2, matches,rtt_dismatches,stc_dismatches, int(sys.argv[3]), _type)
    py.multi_ploter(rtt_name, protein2, matches, int(sys.argv[3]), folder)

#make graphs ------------------------------------------------------------------------------------------------------------------------
    folder = OSfunct.create_dir("../Graphs",rtt_name,protein2,_type)

    if 'd' in _type:
        vtr.freq_VMD(matches,int(sys.argv[3]),sys.argv[8],folder)
    else:
        vtr.freq_VMD(matches,int(sys.argv[3]),"x",folder)

#write recent parameters ------------------------------------------------------------------------------------------------------------
    out.write(protein1[protein1.rfind("/")+1:]+"\n")
    if "/" != sys.argv[4]:
        out.write(sys.argv[4]+"\n")
    else:
        out.write("All\n")
    if "/" != sys.argv[5]:
        out.write(sys.argv[5]+"\n")
    else:
        out.write("All\n")
    out.write(protein2[protein2.rfind("/")+1:]+"\n")
    if "/" != sys.argv[6]:
        out.write(sys.argv[6]+"\n")
    else:
        out.write("All\n")
    if "/" != sys.argv[7]:
        out.write(sys.argv[7]+"\n")
    else:
        out.write("All\n")
    out.write(sys.argv[3]+"\n")

    if "d" in _type:
        out.write("Detailed-")
    else:
        out.write("Simple  -")

    out.write("Distance\n")
    out.close()

    output = open(outname, 'r')
    lines = output.readlines()
    for data in range(0,4):
        print(lines[data])

    
main()
