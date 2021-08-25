import glob
import sys
import time
import re
import argparse
import os
import numpy

cwd = os.getcwd()

parser = argparse.ArgumentParser()
parser.add_argument("identifier", help="number of the tetrazine SM+TS pair to be calculated")
args = parser.parse_args()
tz_number=int(args.identifier)


keepandflag= False # reports if there were maybe problems, does not delete files and keeps them for manual review

#check if confsearch. Use reactand

input=open("xyz/T{0:04d}.xyz".format(tz_number), "r")
lines=input.readlines()
if lines[1].split()[1] == "1":
    confsearch_flag = True
else:
    confsearch_flag = False
    
#check charge
charge=lines[1].split()[0] 
    
input.close()

if confsearch_flag:
#confsearch reactant
    os.system("mkdir T{0:04d}_conf".format(tz_number))
    os.system("cp xyz/T{0:04d}.xyz T{0:04d}_conf/".format(tz_number,tz_number))
    os.chdir("T{0:04d}_conf".format(tz_number))
    # add constraints file because these tight binding methods think planar Tz is ~2 kcal/mol higher than slight tub conformation
    # would double (or triple) the amount of conformers. To prevent this freeze core dihedrals
    # makes conf search slower but QM calcs faster
    output=open(".constrains", "w")
    output.write("$constrain\ndihedral: 3,4,2,5, 0\ndihedral: 4,2,5,6, 0\ndihedral: 5,6,1,3, 0\n force constant = 10\n$end")
    output.close()
    os.system("crest T{0:04d}.xyz -T 20 --chrg {1} > crest.log".format(tz_number,charge))
    log=open("crest.log")
    if "CREST terminated normally" in log.read():
        normalterm= True
    else:
        normalterm=False
    os.chdir("..")    
    #EXIT IF NO NORMALTERM

#confsearch TS
    os.system("mkdir TS_T{0:04d}_conf".format(tz_number))
    os.system("cp xyz/TS_T{0:04d}.xyz TS_T{0:04d}_conf/".format(tz_number,tz_number))
    os.chdir("TS_T{0:04d}_conf".format(tz_number))
    # freeze forming bond length
    # also freeze variables in the TCo to not access higher energy TCO conformers (would be x3)
    output=open(".constrains", "w")
    output.write("$constrain\ndistance: 2,3, 2.24542\ndistance: 1,5,2.24541\ndihedral: 5,11,17,26, 50.992\ndihedral: 17,26,23,20, 117.604\ndihedral: 14,20,23,26, -82.804\nforce constant = 10\n$end")
    output.close()
    os.system("crest TS_T{0:04d}.xyz -T 20 > crest.log".format(tz_number))
    log=open("crest.log")
    if "CREST terminated normally" in log.read():
        normalterm= True
    else:
        normalterm=False
    os.chdir("..")  
    # EXIT IF NO NORMALTERM

    # Now split up xyz in separate files, then convert to g16, then run in g16
    
    # for tetrazine
    os.chdir("T{0:04d}_conf".format(tz_number))
    os.system("mkdir g16")
    os.system("cp crest_conformers.xyz g16/")
    os.chdir("g16")
    input=open("crest_conformers.xyz", "r") #opens multi xyz
    lines=input.readlines() #reads input file and makes it into list
    natoms=int(lines[0]) #finds numebr of atoms for each structure
    nsteps=len(lines)/(natoms+2) #finds number of multiple structures
    input.close()
    for i in range(int(nsteps)):#copy each step in it's own file starts with 0!
        outputname=("T{0:04d}_conf_{1}.xyz".format(tz_number,i))
        output=open(outputname, "w")
        for element in lines[i*(natoms+2):(i+1)*(natoms+2)]:
            output.write(element)
        output.close()
        os.system("python3 ../../G_from_xyz.py T{0:04d}_conf_{1}.xyz {2}".format(tz_number,i,charge))
    os.system("rm crest_conformers.xyz")
    for i in range(int(nsteps)):#copy each step in it's own file starts with 0!
        os.system("g16 T{0:04d}_conf_{1}.com".format(tz_number,i))
        log=open("T{0:04d}_conf_{1}.log".format(tz_number,i))
        if "Normal termination of Gaussian 16" in log.readlines()[-1]:
            pass
        else:
            keepandflag= True
        log.close()
        #check negative frequencies
        log=open("T{0:04d}_conf_{1}.log".format(tz_number,i))
        if "imaginary frequencies" in log.read():
            keepandflag= True
        log.close()
    os.system("python3 ../../GoodVibes_sTRG.py -qh truhlar *.log") # uses goodvibes with slight modification in output
    os.chdir("../..")  
    
    #for TS
    
    os.chdir("TS_T{0:04d}_conf".format(tz_number))
    os.system("mkdir g16")
    os.system("cp crest_conformers.xyz g16/")
    os.chdir("g16")
    input=open("crest_conformers.xyz", "r") #opens multi xyz
    lines=input.readlines() #reads input file and makes it into list
    natoms=int(lines[0]) #finds numebr of atoms for each structure
    nsteps=len(lines)/(natoms+2) #finds number of multiple structures
    input.close()
    for i in range(int(nsteps)):#copy each step in it's own file starts with 0!
        outputname=("TS_T{0:04d}_conf_{1}.xyz".format(tz_number,i))
        output=open(outputname, "w")
        for element in lines[i*(natoms+2):(i+1)*(natoms+2)]:
            output.write(element)
        output.close()
        os.system("python3 ../../G_TS_from_xyz.py TS_T{0:04d}_conf_{1}.xyz {2}".format(tz_number,i,charge))
    os.system("rm crest_conformers.xyz")
    for i in range(int(nsteps)):#copy each step in it's own file starts with 0!
        os.system("g16 TS_T{0:04d}_conf_{1}.com".format(tz_number,i))
        log=open("TS_T{0:04d}_conf_{1}.log".format(tz_number,i))
        if "Normal termination of Gaussian 16" in log.readlines()[-1]:
            pass
        else:
            keepandflag= True
        log.close()
        #check negative frequencies
        log=open("TS_T{0:04d}_conf_{1}.log".format(tz_number,i))
        if "1 imaginary frequencies" in log.read():
            pass
        else:
            keepandflag= True
        log.close()
    os.system("python3 ../../GoodVibes_sTRG.py -qh truhlar *.log") # uses goodvibes with slight modification in output
    os.chdir("../..")
    
   
#now do the same without conf search, use same naming scheme
else:
    #for Tz
    os.system("mkdir T{0:04d}_conf".format(tz_number))
    os.system("mkdir T{0:04d}_conf/g16".format(tz_number))
    os.system("cp xyz/T{0:04d}.xyz T{0:04d}_conf/g16/".format(tz_number,tz_number))
    os.chdir("T{0:04d}_conf/g16/".format(tz_number))
    os.system("mv T{0:04d}.xyz T{0:04d}_conf_0.xyz".format(tz_number))
    os.system("python3 ../../G_from_xyz.py T{0:04d}_conf_0.xyz {1}".format(tz_number,charge))
    os.system("g16 T{0:04d}_conf_0.com".format(tz_number))
    log=open("T{0:04d}_conf_0.log".format(tz_number))
    if "Normal termination of Gaussian 16" in log.readlines()[-1]:
        pass
    else:
        keepandflag= True
    log.close()
    #check negative frequencies
    log=open("T{0:04d}_conf_0.log".format(tz_number))
    if "imaginary frequencies" in log.read():
        keepandflag= True
    log.close()
    os.system("python3 ../../GoodVibes_sTRG.py -qh truhlar *.log") # uses goodvibes with slight modification in output
    os.chdir("../..")
    
    # for TS
    os.system("mkdir TS_T{0:04d}_conf".format(tz_number))
    os.system("mkdir TS_T{0:04d}_conf/g16".format(tz_number))
    os.system("cp xyz/TS_T{0:04d}.xyz TS_T{0:04d}_conf/g16/".format(tz_number,tz_number))
    os.chdir("TS_T{0:04d}_conf/g16/".format(tz_number))
    os.system("mv TS_T{0:04d}.xyz TS_T{0:04d}_conf_0.xyz".format(tz_number))
    os.system("python3 ../../G_TS_from_xyz.py TS_T{0:04d}_conf_0.xyz {1}".format(tz_number,charge))
    os.system("g16 TS_T{0:04d}_conf_0.com".format(tz_number))
    log=open("TS_T{0:04d}_conf_0.log".format(tz_number))
    if "Normal termination of Gaussian 16" in log.readlines()[-1]:
        pass
    else:
        keepandflag= True
    log.close()
    #check negative frequencies
    log=open("TS_T{0:04d}_conf_0.log".format(tz_number))
    if "1 imaginary frequencies" in log.read():
        pass
    else:
        keepandflag= True
    log.close()
    os.system("python3 ../../GoodVibes_sTRG.py -qh truhlar *.log") # uses goodvibes with slight modification in output
    os.chdir("../..")    
    
# now looking for lowest energy conformer, copy them in results folder, cleanup unless review needed

if keepandflag:
    os.system("touch T{0:04d}_review".format(tz_number))
else:
    os.system("mkdir results")
    gibbs = numpy.genfromtxt("T{0:04d}_conf/g16/Goodvibes_output.dat".format(tz_number), usecols=(-1)) # gets qh corrected gibbs free energies
    conf_names = numpy.genfromtxt("T{0:04d}_conf/g16/Goodvibes_output.dat".format(tz_number), dtype=None, encoding=None, usecols=(0)) #gets names
    if gibbs.ndim == 0:
        minimum_e = str(conf_names)
    else:
        minimum_e = conf_names[numpy.where(gibbs == numpy.amin(gibbs))[0]][0]
    os.system("cp T{0:04d}_conf/g16/{1}.log results/T{0:04d}.log".format(tz_number,minimum_e))
        #REMOVE EVERYTHING!
    os.system("rm -r T{0:04d}_conf".format(tz_number))
    
    os.system("mkdir results")
    gibbs = numpy.genfromtxt("TS_T{0:04d}_conf/g16/Goodvibes_output.dat".format(tz_number), usecols=(-1)) # gets qh corrected gibbs free energies
    conf_names = numpy.genfromtxt("TS_T{0:04d}_conf/g16/Goodvibes_output.dat".format(tz_number), dtype=None, encoding=None, usecols=(0)) #gets names
    if gibbs.ndim == 0:
        minimum_e = str(conf_names)
    else:
        minimum_e = conf_names[numpy.where(gibbs == numpy.amin(gibbs))[0]][0]
    os.system("cp TS_T{0:04d}_conf/g16/{1}.log results/TS_T{0:04d}.log".format(tz_number,minimum_e))
        #REMOVE EVERYTHING!
    os.system("rm -r TS_T{0:04d}_conf".format(tz_number))

