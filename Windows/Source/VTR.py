# VTR 2.0 for Windows
import os

os.system("md ..\\Results")
os.system("md ..\\Results\\Contacts")
os.system("md ..\\Results\\Matches")
os.system("md ..\\Results\\Dismatches")
os.system("md ..\\Plots")
os.system("md ..\\Graphs")
while (1):
     os.system("cls")
     os.system("echo ____________________________________________________________\n")
     os.system("echo VTR Geometric v 2.0\n")
     os.system("echo ____________________________________________________________\n")
     os.system("echo Proteins in directory: ")
     os.system("dir ..\Data /a-d /on")
     os.system("echo  ___________________________________________________________\n")
     rtt_protein = input("Choose the first:")
     while (not(os.path.exists("..\\Data\\" + rtt_protein + ".pdb"))):
          print("This file does not exist...")
          rtt_protein = input("Choose the first:")
     stc_protein = input("Choose the second:")
     while (not(os.path.exists("..\\Data\\" + stc_protein + ".pdb"))):
          print("This file does not exist...")
          stc_protein = input("Choose the second:")
     cutoff = input("Choose the cutoff:")
     chainrtt1 = input("Choose the first chain for " + rtt_protein +"(A, B ..., press enter for run all x all chains):")
     if chainrtt1 == "":
          chainrtt1 = "/"
     chainrtt2 = input("Choose the second chain for " + rtt_protein +"(A, B ..., press enter for run all x all chains):")
     if chainrtt2 == "":
          chainrtt2 = "/"
     chainstc1 = input("Choose the first chain for " + stc_protein +"(A, B ..., press enter for run all x all chains):")
     if chainstc1 == "":
          chainstc1 = "/"
     chainstc2 = input("Choose the second chain for " + stc_protein +"(A, B ..., press enter for run all x all chains):")
     if chainstc2 == "":
          chainstc2 = "/"

     question = 0
     selection = '''Please Select the Execution Type:
     (1) Simple Graphics, Match by Distance
     (2) Detailed Graphics, Match by Distance'''
     
     while (question > 2 or question < 1):
          print(selection)
          question = int(input('   '))
     if (question == 1):
          execute = "python VTR_Geometric.py ../Data/" + rtt_protein + ".pdb ../Data/" + stc_protein + ".pdb " + cutoff + " " + chainrtt1 + " " + chainrtt2 + " " + chainstc1 + " " + chainstc2
     elif (question == 2):
          execute = "python VTR_Geometric.py ../Data/" + rtt_protein + ".pdb ../Data/" + stc_protein + ".pdb " + cutoff + " " + chainrtt1 + " " + chainrtt2 + " " + chainstc1 + " " + chainstc2 + " d"
     os.system(execute)
     input('Press Enter to Return')
     

