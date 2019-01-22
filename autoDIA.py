#!/usr/bin/env python 

"""
*****************************************************************************************
                                    autoDIA
								Dennis Svatunek
								    UCLA
							d.svatunek@chem.ucla.com
							
						Download the newest version from:
					   https://github.com/dsvatunek/autoDIA
					   
					   Please cite:
*****************************************************************************************				   

Python wrapper for Gaussian to automatically perform Distortion/Interaction analysis calculations (also know as Activation Strain Analysis).

Further reading:

"""

__version__= ' '

#Python2 compatibility
try:
	range = xrange
except NameError:
	pass

radian2degree= 57.2958
hartree2kcal = 627.509

#Periodic table
periodic_table = ["","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
    "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl",
    "Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"]

# checks if "s" is an int, returns true or false
def isInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False	

def main():
	#getting input file
	parser =  argparse.ArgumentParser(usage='%(prog)s inp_file',  description=description, epilog=epilog)
	parser.add_argument('inp_file', metavar='Input file', type=str, help='Name of the input file')
	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit(1)	
	args = parser.parse_args()
	# open input file
	input_object = open(args.inp_file, 'r')
#----------------------------------------------------------------------------------------
	overallstarttime=time.time()	
#-------------------------------------Parse Settings and Inputs----------------------
	starttime=time.time() # start timing for input parsing
	settings, structures = parse_in()
	endtime=time.time()
	totaltime=str(endtime-starttime)
	seconds=totaltime.split('.')[0]
	milliseconds=float('0.'+totaltime.split('.')[1])*1000
	log = open(logfile, 'a')
	log.write('\nProcessed input in {} seconds and {:.0f} ms\n'.format(seconds, float(milliseconds)))
	log.write('-------------------------------------------------------\n\n')
	log.close()
#----------------------------------------------------------------------------------------
	
	
	
	

#-------------------------------------Clean up-----------------------------------------
	
		


if __name__ == "__main__":
    main()