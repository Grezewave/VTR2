#This program recieve a PDB format file and create a protein object
class proteins:
    def __init__(self):
        self.idPDB = ""
        self.header = ""
        self.title = ""
        self.chains = []
    def size(self):
        size = 0
        for i in self.chains:
            for e in i.residues:
                size += len(e.atoms)
        return size
class chain:
    def __init__(self):
        self.id = ""
        self.residues = []
class residue:
    def __init__(self):
        self.id = ""
        self.parameter = 0
        self.atoms = []
class atom:
    def __init__(self):
        self.id = ""
        self.type = ""
        self.x = 0
        self.y = 0
        self.z = 0
        self.occupancy = 0
        self.b_factor = 0
        
def tellme(text,data):
    #Return the line with the "text" identifier
    while(len(text)<6):
        text += " "
    j = 6
    info = ""
    #search for two spaces in sequence, to stop the information collection
    for i in data:
        if i[0:6] == text:
            while i[j] == " ":
                j+=1
            k = j
            while i[k+1] != "\n":
                k+=1
            info += i[j:k]
    return info

def idPDB(data):
    #Return the protein ID
    return data[0][62:66]
        
def chainslist(data):
    #Return a list with all chains based in ATOM section
    info = []
    for i in data:
        if i[0:4] == "ATOM":
            if not(i[21] in info):
                info.extend(i[21])
    return info
            
def chainsdef(chainlist,data):
    #Create a list of "chain" objects
    info = ""
    e = 0
    guard = []
    #Turn each element of the chains list in a object, and put id
    for i in range(0,len(chainlist)):
        info = chainlist[i]
        chainlist[i] = chain()
        chainlist[i].id = info
    #Create the residues list of each chain
    for i in data:
        if i[0:4] == "ATOM":
            if i[21] == chainlist[e].id:
                if int(i[22:26]) not in guard:
                    guard.append(int(i[22:26]))
                    chainlist[e].residues.append(i[17:20])
            else:
                e+=1
                guard = []
                if i[21] == chainlist[e].id:
                    if int(i[22:26]) not in guard:
                        guard.append(int(i[22:26]))
                        chainlist[e].residues.append(i[17:20])
    return chainlist

def residuedef(data,reslist,chain):
    #Create a list of "residues" objects
    o = -1
    info = reslist[:]
    guard = []
    #Turn each element of the residues list in a object, and put id
    for i in range(0,len(reslist)):
        reslist[i] = residue()
        reslist[i].id = info[i]
    #Put the other informations of the each residue, and create a list of atoms for eache residue
    for i in data:
        if i[0:4] == "ATOM":
            if i[21] == chain:
                if int(i[22:26]) not in guard:
                    guard.append(int(i[22:26]))
                    o+=1
                    if i[17:20] == reslist[o].id:
                        reslist[o].parameter = int(i[22:26])
                        reslist[o].atoms.append(i[13:17])
                else:
                    if i[17:20] == reslist[o].id:
                        reslist[o].parameter = int(i[22:26])
                        reslist[o].atoms.append(i[13:17])
    return reslist
    
def proteindef(data):
    #Create the protein object
    protein = proteins()
    protein.idPDB = idPDB(data)
    protein.header = tellme("HEADER",data)
    protein.title = tellme("TITLE",data)
    protein.chains = chainslist(data)
    protein.chains = chainsdef(protein.chains,data)
    for i in protein.chains:
        i.seq = residuedef(data,i.residues,i.id)
        for e in i.residues:
            e.atoms = atomdef(data,e.atoms,e.id,e.parameter,i.id)
    return protein

                
def atomdef(data,atomlist,resname,parameter,chain):
    #Create a list of "atoms" objects
    info = atomlist[:]
    x = 0
    #Turn each element of the atoms list in a object, and put id
    for e in range(0,len(atomlist)):
        atomlist[e] = atom()
        atomlist[e].type = info[e]
    #Put the other informations of the each atom
    for i in data:
        if i[0:4] == "ATOM":
            if i[17:20] == resname and i[21] == chain and int(i[22:26]) == parameter:
                if i[13:17] == atomlist[x].type:
                    atomlist[x].id = int(i[6:11])
                    atomlist[x].x = float(i[30:38])
                    atomlist[x].y = float(i[38:46])
                    atomlist[x].z = float(i[46:54])
                    try:
                        atomlist[x].occupancy = float(i[54:60])
                    except:
                        h = 0
                    try:
                        atomlist[x].b_factor = float(i[60:66])
                    except:
                        h = 0
                    x+=1
    return atomlist

def classify(file):
    #Read file and create object
    reader = open(file,'r')
    data = reader.readlines()
    reader.close()
    protein = proteindef(data)
    protein.title = protein.title.replace("  ","")
    return protein
                
    
