import glob
import sys
import time
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("filename", help="xyz filename")
parser.add_argument("charge", help="xyz charge")
args = parser.parse_args()
file=args.filename
charge=args.charge


def createcom(name): #makes *.txt file with input name, exits py if not possible

        try:    #creating name.txt
                output = open(name + '.com','w')   # Trying to create a new file or open one
                output.close()
        except:
                print('Something went wrong! Couldn\'t create output.txt!')
                time.sleep(5) # to see what happened
                sys.exit(0) # quit Python


input=open(file, "r")
createcom(file[0:-4])
output=open(file[0:-4]+".com", "a")
output.write("%mem=56GB\n%nprocshared=20\n# opt freq wb97xd def2SVP scrf=(smd)\n\n[No Title]\n\n")
output.write("{0} 1\n".format(charge))
lines=input.readlines() 
lines=lines[2:]
for i in range(len(lines)):
    output.write(lines[int(i)])
output.write("\n")
output.close()