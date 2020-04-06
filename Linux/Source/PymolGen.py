import OSfunct
import os
import Plot

def colorchange(index):
    #return the color of pymol type
    colors = ["aquamarine","blue","bluewhite","br0","br1","br2","br3","br4","br5","br6","br7","br8","br9","brightorange","brown","carbon","chartreuse","chocolate","cyan","darksalmon","dash","deepblue","deepolive","deeppurple","deepsalmon","deepteal","density","dirtyviolet","firebrick","forest","gray","green","greencyan","grey","hotpink","hydrogen","lightblue","lightmagenta","lightorange","lightpink","lightteal","lime","limegreen","limon","magenta","marine","nitrogen","olive","orange","oxygen","palecyan","palegreen","paleyellow","pink","purple","purpleblue","raspberry","red","ruby","salmon","sand","skyblue","slate","smudge","splitpea","sulfur","teal","tv_blue","tv_green","tv_orange","tv_red","tv_yellow","violet","violetpurple","warmpink","wheat","white","yellow","yelloworange"]
    return (colors[(index)%(len(colors)-1)])

def detailed_ploter(rtt_path, stc_path, matches, rtt_dismatches, stc_dismatches, cutoff, _type):
    #create a pml file for the full structure superposition with the matchs in highlights
    folder = OSfunct.create_dir("../Plots",rtt_path,stc_path,_type)
    #define ne file name
    pmlname = "../Plots/"+folder+"/c_scale" + rtt_path[rtt_path.rfind("/")+1:rtt_path.rfind("_")] + "_x_" + stc_path[stc_path.rfind("/")+1:-4] + ".pml"
    pml = open(pmlname,'w')
    pml.write("load " + rtt_path[rtt_path.rfind("/")+1:] + "\n")
    pml.write("load " + stc_path[stc_path.rfind("/")+1:] + "\n")
    x = 0
    #put each match in the visualization, with a color gradded equivalent to VMD
    for i in matches:
        color = Plot.colorscale(i.VMD(),cutoff,'l')
        entry = "set_color " + str(color[0]) + "_" + str(color[1]) + "_" + str(color[2]) + ", [" + str(color[0]) + "," + str(color[1]) + "," + str(color[2]) + "]\n"
        pml.write(entry)
        selection1 = rtt_path[rtt_path.rfind("/")+1:rtt_path.rfind("_")] + str(i.rtt_contact.atom1.id)
        entry = "select " + selection1 + ",model " + rtt_path[rtt_path.rfind("/")+1:-4] + " and id " + str(i.rtt_contact.atom1.id) + "\n"
        pml.write(entry)
        selection2 = rtt_path[rtt_path.rfind("/")+1:rtt_path.rfind("_")] + str(i.rtt_contact.atom2.id)
        entry = "select " + selection2 + ",model " + rtt_path[rtt_path.rfind("/")+1:-4] + " and id " + str(i.rtt_contact.atom2.id) + "\n"
        pml.write(entry)
        entry = "distance " + selection1 + "-" + selection2 + ", " + selection1 + ", " + selection2 + "\n"
        pml.write(entry)
        entry = "color " + str(color[0]) + "_" + str(color[1]) + "_" + str(color[2]) + ", " + selection1 + "-" + selection2 + "\n"
        pml.write(entry)
        selection1 = stc_path[stc_path.rfind("/")+1:-4] + str(i.stc_contact.atom1.id)
        entry = "select " + selection1 + ",model " + stc_path[stc_path.rfind("/")+1:-4] + " and id " + str(i.stc_contact.atom1.id) + "\n"
        pml.write(entry)
        selection2 = stc_path[stc_path.rfind("/")+1:-4] + str(i.stc_contact.atom2.id)
        entry = "select " + selection2 + ",model " + stc_path[stc_path.rfind("/")+1:-4] + " and id " + str(i.stc_contact.atom2.id) + "\n"
        pml.write(entry)
        entry = "distance " + selection1 + "-" + selection2 + ", " + selection1 + ", " + selection2 + "\n"
        pml.write(entry)
        entry = "color " + str(color[0]) + "_" + str(color[1]) + "_" + str(color[2]) + ", " + selection1 + "-" + selection2 + "\n"
        pml.write(entry)
        x += 1
    pml.write("set_color 125_125_125,[125,125,125]\n")
    entry = ' and('
    #Put the dismatches in a separated sellection,on gray color
    for i in rtt_dismatches:
        entry+=" id "+str(i.atom1.id)
        entry+=" or id "+str(i.atom2.id)+" or "
    pml.write('select '+rtt_path[rtt_path.rfind("/")+1:rtt_path.rfind("_")]+'_dismatches,model '+rtt_path[rtt_path.rfind("/")+1:rtt_path.rfind("_")]+entry[:-4]+')\n')
    pml.write("color 125_125_125, "+rtt_path[rtt_path.rfind("/")+1:rtt_path.rfind("_")]+"_dismatches\n")

    entry = ' and('
    for i in stc_dismatches:
        entry+=" id "+str(i.atom1.id)
        entry+=" or id "+str(i.atom2.id)+" or "
    pml.write('select '+stc_path[stc_path.rfind("/")+1:-4]+'_dismatches,model '+stc_path[stc_path.rfind("/")+1:-4]+entry[:-4]+')\n')
    pml.write("color 125_125_125, "+stc_path[stc_path.rfind("/")+1:-4]+"_dismatches\n")
    
    hideH = "sele resn HOH\nhide (sele)" 
    pml.write(hideH)
    pml.close()
    #return the folder name
    return(folder)

def multi_ploter(rtt_path, stc_path, matches, cutoff, folder):
    #create a visualization of each contact on Pymol
    for i in matches:
        #load the proteins, and hide everything
        pmlname = "../Plots/" + folder + "/" + i.rtt_contact.residue1.id + str(i.rtt_contact.residue1.parameter) + "--" + i.rtt_contact.residue2.id + str(i.rtt_contact.residue2.parameter) + "_x_" + i.stc_contact.residue1.id + str(i.stc_contact.residue1.parameter) + "--" + i.stc_contact.residue2.id + str(i.stc_contact.residue2.parameter) + ".pml"
        pml = open(pmlname,'a')
        pml.write("load " + rtt_path[rtt_path.rfind("/")+1:] + "\n")
        pml.write("load " + stc_path[stc_path.rfind("/")+1:] + "\n")
        entry = "hide all\n"
        pml.write(entry)
        #show residues as sticks
        selection1 = rtt_path[rtt_path.rfind("/")+1:rtt_path.rfind("_")] + i.rtt_contact.residue1.id + str(i.rtt_contact.residue1.parameter)
        entry = "select " + selection1 + ",model " + rtt_path[rtt_path.rfind("/")+1:-4] + " and chain " + i.rtt_contact.chain1.id + " and resi " + str(i.rtt_contact.residue1.parameter) + "\n"
        pml.write(entry)
        selection2 = rtt_path[rtt_path.rfind("/")+1:rtt_path.rfind("_")] + i.rtt_contact.residue2.id + str(i.rtt_contact.residue2.parameter)
        entry = "select " + selection2 + ",model " + rtt_path[rtt_path.rfind("/")+1:-4] + " and chain " + i.rtt_contact.chain2.id + " and resi " + str(i.rtt_contact.residue2.parameter) + "\n"
        pml.write(entry)
        entry = "show sticks, " + selection1 + " " + selection2 + "\n"
        pml.write(entry)
        selection1 = stc_path[stc_path.rfind("/")+1:-4] + i.stc_contact.residue1.id + str(i.stc_contact.residue1.parameter)
        entry = "select " + selection1 + ",model " + stc_path[stc_path.rfind("/")+1:-4] + " and chain " + i.stc_contact.chain1.id + " and resi " + str(i.stc_contact.residue1.parameter) + "\n"
        pml.write(entry)
        selection2 = stc_path[stc_path.rfind("/")+1:-4] + i.stc_contact.residue2.id + str(i.stc_contact.residue2.parameter)
        entry = "select " + selection2 + ",model " + stc_path[stc_path.rfind("/")+1:-4] + " and chain " + i.stc_contact.chain2.id + " and resi " + str(i.stc_contact.residue2.parameter) + "\n"
        pml.write(entry)
        entry = "show sticks, " + selection1 + " " + selection2 + "\n"
        pml.write(entry)
        entry = "sele resn HOH\nhide (sele)\n"
        pml.write(entry)
    for i in matches:
        #open the same file
        pmlname = "../Plots/" + folder + "/" + i.rtt_contact.residue1.id + str(i.rtt_contact.residue1.parameter) + "--" + i.rtt_contact.residue2.id + str(i.rtt_contact.residue2.parameter) + "_x_" + i.stc_contact.residue1.id + str(i.stc_contact.residue1.parameter) + "--" + i.stc_contact.residue2.id + str(i.stc_contact.residue2.parameter) + ".pml"
        pml = open(pmlname,'a')
        #set each match, criating a dashed line between atoms
        color = Plot.colorscale(i.VMD(),cutoff,'l')
        entry = "set_color " + str(color[0]) + "_" + str(color[1]) + "_" + str(color[2]) + ", [" + str(color[0]) + "," + str(color[1]) + "," + str(color[2]) + "]\n"
        pml.write(entry)
        selection1 = rtt_path[rtt_path.rfind("/")+1:rtt_path.rfind("_")] + str(i.rtt_contact.atom1.id)
        entry = "select " + selection1 + ",model " + rtt_path[rtt_path.rfind("/")+1:-4] + " and id " + str(i.rtt_contact.atom1.id) + "\n"
        pml.write(entry)
        selection2 = rtt_path[rtt_path.rfind("/")+1:rtt_path.rfind("_")] + str(i.rtt_contact.atom2.id)
        entry = "select " + selection2 + ",model " + rtt_path[rtt_path.rfind("/")+1:-4] + " and id " + str(i.rtt_contact.atom2.id) + "\n"
        pml.write(entry)
        entry = "distance " + selection1 + "-" + selection2 + ", " + selection1 + ", " + selection2 + "\n"
        pml.write(entry)
        entry = "color " + str(color[0]) + "_" + str(color[1]) + "_" + str(color[2]) + ", " + selection1 + "-" + selection2 + "\n"
        pml.write(entry)
        selection1 = stc_path[stc_path.rfind("/")+1:-4] + str(i.stc_contact.atom1.id)
        entry = "select " + selection1 + ",model " + stc_path[stc_path.rfind("/")+1:-4] + " and id " + str(i.stc_contact.atom1.id) + "\n"
        pml.write(entry)
        selection2 = stc_path[stc_path.rfind("/")+1:-4] + str(i.stc_contact.atom2.id)
        entry = "select " + selection2 + ",model " + stc_path[stc_path.rfind("/")+1:-4] + " and id " + str(i.stc_contact.atom2.id) + "\n"
        pml.write(entry)
        entry = "distance " + selection1 + "-" + selection2 + ", " + selection1 + ", " + selection2 + "\n"
        pml.write(entry)
        entry = "color " + str(color[0]) + "_" + str(color[1]) + "_" + str(color[2]) + ", " + selection1 + "-" + selection2 + "\n"
        pml.write(entry)
