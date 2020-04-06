#Contain the Operational System's functions

import os
import time
from datetime import datetime

def TMAlign(protein1,protein2):
    #Delete existing alignments and create new
    if (not(os.path.exists("tmalign.exe"))):
        os.system("g++ TMAlign.cpp -o tmalign")
        print("TMAlign compilado!")

    path = protein1[protein1.rfind("/")+1:-4] + "x" + protein2[protein2.rfind("/")+1:-4] + "_align"
    if os.path.exists("../Data/" + path):
        os.system("rd -r ..\\Data\\" + path + "/s /q")
    os.system("md ..\\Data\\" + path)
    callalign = "tmalign " + protein1 + " " + protein2 + " -o " + "../Data/" + path + "/" + protein1[protein1.rfind("/")+1:-4]
    os.system(callalign)
    return(path)

def create_dir(location,rtt_path,stc_path,_type):
    #create a directory in location
    folder = rtt_path[rtt_path.rfind("/")+1:rtt_path.rfind("_")] + "_x_" + stc_path[stc_path.rfind("/")+1:-4] + _type + "_plot"
    if (os.path.exists(location + "/" + folder)):
        pmlname = "rd " + location.replace("/","\\") + "\\" + folder + "/s /q"
        os.system(pmlname)
    pmlname = "md " + location.replace("/","\\") + "\\" + folder
    os.system(pmlname)
    return(folder)
