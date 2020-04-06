#VTR Functions

import Contacts
import matplotlib
import matplotlib.pyplot as plt
import numpy
import OSfunct
import math

class match:
    def __init__(self,rtt_contact,stc_contact):
        self.rtt_contact = rtt_contact
        self.stc_contact = stc_contact
    def Vector11(self):
        return (Contacts.adistance(self.rtt_contact.atom1,self.stc_contact.atom1))
    def Vector12(self):
        return (Contacts.adistance(self.rtt_contact.atom1,self.stc_contact.atom2))
    def Vector21(self):
        return (Contacts.adistance(self.rtt_contact.atom2,self.stc_contact.atom1))
    def Vector22(self):
        return (Contacts.adistance(self.rtt_contact.atom2,self.stc_contact.atom2))
    def VMD(self):
        VMD=[]
        VMD.append((self.Vector11() + self.Vector22())/2)
        VMD.append((self.Vector12() + self.Vector21())/2)
        return (min(VMD))
            
            

def minVMD(matches,blacklist,cutoff):
    #Find the lower VMD in the list of possible matches and remove invalid matches from the list
    minVMD = cutoff
    #cria uma cópia da lista
    temp_matches = matches
    match = 0
    for i in matches:
        #search for the lower VMD, if members are not in blacklist
        if not((i.rtt_contact in blacklist) or (i.stc_contact in blacklist)):
            if (i.VMD() <= minVMD):
                match = i
                minVMD = match.VMD()
        else:
            #remove invalid from list, this improve the eficience of futures searchs
            temp_matches.remove(i)
    return (match,temp_matches)
    
def RMSD(matches, protein1, protein2):
    RMSD = 0
    for i in matches:
        RMSD += (i.VMD()**2)
    RMSD = (RMSD/len(matches))**(1/2)
    return RMSD

def VTR(matches, protein1, protein2, rtt_dismatches, stc_dismatches):
    VTR = 0
    for i in matches:
        VTR += i.VMD()
    VTR = ((VTR/len(matches))*(1+abs(protein1.size()-protein2.size()))) + rtt_dismatches + stc_dismatches
    return VTR

def match_contacts_equivalence(rtt_contacts,stc_contacts,cutoff,equivalences):
    #try to match contacts from a list of equivalences, output of Multiprot
    all_matches = []
    rtt_dismatches = []
    stc_dismatches = []
    blacklist = []
    matches = []
    create = 0

    for i in rtt_contacts:
        for j in stc_contacts:
            try:
                #verify if the current residues of contacts are equivalent
                if [i.chain1.id,j.chain1.id] in equivalences[(i.residue1.parameter,j.residue1.parameter)]:
                    if [i.chain2.id,j.chain2.id] in equivalences[(i.residue2.parameter,j.residue2.parameter)]:
                        create = 1
                elif [i.chain1.id,j.chain2.id] in equivalences[(i.residue1.parameter,j.residue2.parameter)]:             
                    if [i.chain2.id,j.chain2.id] in equivalences[(i.residue2.parameter,j.residue2.parameter)]:
                        create = 1
            except:
                r = 0
            #if true, create a possible contact
            if create:
                _match = match(i,j)
                if _match.VMD() <= cutoff:
                    all_matches.append(_match)

            create = 0

    #refine the contacts, serching for the bests    
    control = True
    while control:
        _match,all_matches = minVMD(all_matches, blacklist, cutoff)
        if _match != 0:
            blacklist.append(_match.rtt_contact)
            blacklist.append(_match.stc_contact)
            matches.append(_match)
            all_matches.remove(_match)
        else:
            control = False

    #make the dismatch's lists
    for m in rtt_contacts:
        if not(m in blacklist):
            rtt_dismatches.append(m)
    for m in stc_contacts:
        if not(m in blacklist):
            stc_dismatches.append(m)
    return (matches,rtt_dismatches,stc_dismatches)

def match_contacts(rtt_contacts,stc_contacts,cutoff):
    #try to match contacts by distance only
    all_matches = []
    rtt_dismatches = []
    stc_dismatches = []
    blacklist = []
    matches = []
    for i in rtt_contacts:
        for j in stc_contacts:
            #create a list of possibles contacts
            _match = match(i,j)
            if _match.VMD() <= cutoff:
                all_matches.append(_match)

    #refine the contacts, serching for the bests                
    control = True
    while control:
        _match,all_matches = minVMD(all_matches, blacklist, cutoff)
        if _match != 0:
            blacklist.append(_match.rtt_contact)
            blacklist.append(_match.stc_contact)
            matches.append(_match)
            all_matches.remove(_match)
        else:
            control = False
    
    #make the dismatch's lists
    for m in rtt_contacts:
        if not(m in blacklist):
            rtt_dismatches.append(m)
    for m in stc_contacts:
        if not(m in blacklist):
            stc_dismatches.append(m)
    return (matches,rtt_dismatches,stc_dismatches)

def freq_VMD(matches,cutoff,detail,folder):
    #create histograms for VMD frequency and residues interactions
    x = [i for i in numpy.arange(0,float(cutoff),cutoff/20)]
    y = []
    vmd = []
    for match in matches:
        #create a list of VMDs
        y.append(match.VMD())
    #create lists for histogram 
    frequency, _vmd, thrash = plt.hist(y,bins = x)
    for i in range(1,len(_vmd)):
        #create data for a line graph
        vmd.append((_vmd[i-1] + _vmd[i])/2)

    plt.plot(vmd,frequency)
    plt.title('Frequency distribuition')
    plt.ylabel('Frequency')
    plt.xlabel('VMD')

    plt.savefig('../Graphs/' + folder + '/VMD')
    plt.show()

    if detail == "d":
        #provide detailed graphic info
        #initiate a frequency distribuition for each reasidue in each protein
        rtt_freq = { "ALA" : [0 for i in range(0,20)],
                     "ARG" : [0 for i in range(0,20)],
                     "ASN" : [0 for i in range(0,20)],
                     "ASP" : [0 for i in range(0,20)],
                     "CYS" : [0 for i in range(0,20)],
                     "GLN" : [0 for i in range(0,20)],
                     "GLY" : [0 for i in range(0,20)],
                     "GLU" : [0 for i in range(0,20)],
                     "HIS" : [0 for i in range(0,20)],
                     "ILE" : [0 for i in range(0,20)],
                     "LEU" : [0 for i in range(0,20)],
                     "LYS" : [0 for i in range(0,20)],
                     "MET" : [0 for i in range(0,20)],
                     "PHE" : [0 for i in range(0,20)],
                     "PRO" : [0 for i in range(0,20)],
                     "SER" : [0 for i in range(0,20)],
                     "THR" : [0 for i in range(0,20)],
                     "TRP" : [0 for i in range(0,20)],
                     "TYR" : [0 for i in range(0,20)],
                     "VAL" : [0 for i in range(0,20)]}
                    
        stc_freq = { "ALA" : [0 for i in range(0,20)],
                     "ARG" : [0 for i in range(0,20)],
                     "ASN" : [0 for i in range(0,20)],
                     "ASP" : [0 for i in range(0,20)],
                     "CYS" : [0 for i in range(0,20)],
                     "GLN" : [0 for i in range(0,20)],
                     "GLY" : [0 for i in range(0,20)],
                     "GLU" : [0 for i in range(0,20)],
                     "HIS" : [0 for i in range(0,20)],
                     "ILE" : [0 for i in range(0,20)],
                     "LEU" : [0 for i in range(0,20)],
                     "LYS" : [0 for i in range(0,20)],
                     "MET" : [0 for i in range(0,20)],
                     "PHE" : [0 for i in range(0,20)],
                     "PRO" : [0 for i in range(0,20)],
                     "SER" : [0 for i in range(0,20)],
                     "THR" : [0 for i in range(0,20)],
                     "TRP" : [0 for i in range(0,20)],
                     "TYR" : [0 for i in range(0,20)],
                     "VAL" : [0 for i in range(0,20)]}

        residues = ["ALA","ARG","ASN","ASP","CYS","GLN","GLY","GLU","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]

        x = 1
        #walks the dictionary and the list, counting the frequency of interactions between residues
        for ref in rtt_freq:
            for i in range(0,len(residues)):
                for match in matches:
                    if (match.rtt_contact.residue1.id == ref):
                        if (match.rtt_contact.residue2.id == residues[i]):
                            rtt_freq[ref][i] += 1
                    elif (match.rtt_contact.residue2.id == ref):
                        if (match.rtt_contact.residue1.id == residues[i]):
                            rtt_freq[ref][i] += 1
                    if (match.stc_contact.residue1.id == ref):
                        if (match.stc_contact.residue2.id == residues[i]):
                            stc_freq[ref][i] += 1
                    elif (match.stc_contact.residue2.id == ref):
                        if (match.stc_contact.residue1.id == residues[i]):
                            stc_freq[ref][i] += 1
        _max = 0
        #wals the frequency lists, seraching for the max, for normalization of graphs
        for ref in rtt_freq:
            for index in rtt_freq[ref]:
                if index > _max:
                    _max = index

            for index in stc_freq[ref]:
                if index > _max:
                    _max = index
        #plot graphs
        pic = plt.figure(figsize=(18, 5))
        for ref in rtt_freq:
            plt.subplot(121)
            plt.bar(residues, rtt_freq[ref], width = 0.7)
            plt.axis([-1,20,0,_max])
            plt.subplot(122)
            plt.bar(residues, stc_freq[ref], width = 0.7)
            plt.axis([-1,20,0,_max])
            plt.subplots_adjust(left=0.04, bottom=0.05, right=0.99, top=0.97, wspace=0.09, hspace=0.41)
            plt.savefig('../Graphs/' + folder + '/' + ref)
            plt.clf()
            x += 1
        print('Files can be found in the Graphs folder')

def write_dismatch(protein1,protein2,rtt_dismatches,stc_dismatches,_type):
    #write the dismatche file
    dis = open('../Results/Dismatches/'+protein1[protein1.rfind("/")+1:-4]+"x"+protein2[protein2.rfind("/")+1:-4]+_type+'.txt','w')
    dis.write(protein1[protein1.rfind("/")+1:-4]+' No matched Contacts:\n\n')
    for i in rtt_dismatches:
        dis.write(str(i.number)+' - '+i.chain1.id+" "+i.residue1.id+" "+str(i.residue1.parameter)+" "+i.atom1.type+" VS "+i.chain2.id+" "+i.residue2.id+" "+str(i.residue2.parameter)+" "+i.atom2.type+"\n")
        dis.write("Distance: "+str(i.distance())+" - ")
        for e in i.contacts:
            dis.write(e+", ")
        dis.write("\n")
    
    dis.write("\n\n")
    dis.write(protein2[protein2.rfind("/")+1:-4]+' No matched Contacts:\n\n')
    for i in stc_dismatches:
        dis.write(str(i.number)+' - '+i.chain1.id+" "+i.residue1.id+" "+str(i.residue1.parameter)+" "+i.atom1.type+" VS "+i.chain2.id+" "+i.residue2.id+" "+str(i.residue2.parameter)+" "+i.atom2.type+"\n")
        dis.write("Distance: "+str(i.distance())+" - ")
        for e in i.contacts:
            dis.write(e+", ")
        dis.write("\n")
 
def writer(protein1,protein2,rtt_protein,stc_protein,rtt_contacts,stc_contacts,matches,_type):
    #write the match file
    outfile = "../Results/Matches/" + protein1[protein1.rfind("/")+1:-4] + "x" + protein2[protein2.rfind("/")+1:-4] + _type + ".txt"
    out = open(outfile,'w')
    out.write("Rotate Protein: ")
    out.write(rtt_protein.idPDB)
    out.write("\n")
    out.write(rtt_protein.header)
    out.write("\n\nStatic Protein: ")
    out.write(stc_protein.idPDB)
    out.write("\n")
    out.write(stc_protein.header)
    out.write("\n")
    for i in matches:
        out.write(i.rtt_contact.chain1.id)
        out.write(" ")
        out.write(i.rtt_contact.residue1.id)
        out.write(" ")
        out.write(str(i.rtt_contact.residue1.parameter))
        out.write(" ")
        out.write(i.rtt_contact.atom1.type)
        out.write(" ")
        out.write(str(i.rtt_contact.atom1.id))
        out.write(" ----------- ")
        out.write(i.rtt_contact.chain2.id)
        out.write(" ")
        out.write(i.rtt_contact.residue2.id)
        out.write(" ")
        out.write(str(i.rtt_contact.residue2.parameter))
        out.write(" ")
        out.write(i.rtt_contact.atom2.type)
        out.write(" ")
        out.write(str(i.rtt_contact.atom2.id))
        out.write("\n")
        out.write(i.stc_contact.chain1.id)
        out.write(" ")
        out.write(i.stc_contact.residue1.id)
        out.write(" ")
        out.write(str(i.stc_contact.residue1.parameter))
        out.write(" ")
        out.write(i.stc_contact.atom1.type)
        out.write(" ")
        out.write(str(i.stc_contact.atom1.id))
        out.write(" ----------- ")
        out.write(i.stc_contact.chain2.id)
        out.write(" ")
        out.write(i.stc_contact.residue2.id)
        out.write(" ")
        out.write(str(i.stc_contact.residue2.parameter))
        out.write(" ")
        out.write(i.stc_contact.atom2.type)
        out.write(" ")
        out.write(str(i.stc_contact.atom2.id))
        out.write("\n")
        out.write("VMD: ")
        out.write(str(i.VMD()))
        out.write("            Contact types: ")
        for e in i.rtt_contact.contacts:
            out.write(e)
            out.write(", ")
        out.write(" / ")
        for e in i.stc_contact.contacts:
            out.write(e)
            out.write(", ")    
        out.write("\n\n")
    out.close()
