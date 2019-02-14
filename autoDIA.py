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

from DIA_inputparser import *
import argparse, sys, time, subprocess

__version__= ' '

#Python2 compatibility
try:
	range = xrange
except NameError:
	pass

radian2degree= 57.2958
hartree2kcal = 627.509


def main():
	#getting input file
	parser =  argparse.ArgumentParser(usage='%(prog)s inp_file')
	parser.add_argument('inp_file', metavar='Input file', type=str, help='Name of the input file')
	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit(1)	
	args = parser.parse_args()
#----------------------------------------------------------------------------------------
	overallstarttime=time.time()	
#-------------------------------------Parse Settings and Inputs----------------------
	starttime=time.time() # start timing for input parsing
	print("test")
	settings, structures = parse_in(args.inp_file)
	print("test")
	print(structures.xyz)
	print(dir(settings))
	endtime=time.time()
	totaltime=str(endtime-starttime)
	seconds=totaltime.split('.')[0]
	milliseconds=float('0.'+totaltime.split('.')[1])*1000
	
	log = open(settings.logfile, 'w')
	hostname = subprocess.check_output(['hostname']).decode('utf-8')
	username = subprocess.check_output(['whoami']).decode('utf-8')
	log.write(str(username).strip()+'@'+str(hostname).strip()+'\n')
	log.write(time.ctime()+'\n')
	log.write('-------------------------------------------------------\n')
	log.write('Input file: '+args.inp_file)
	log.write('\nProcessed input and parsed structures in {} seconds and {:.0f} ms\n'.format(seconds, float(milliseconds)))
	log.write('-------------------------------------------------------\n')
	log.write('Determined settings and variables: \n\n')
	log.write('{0:50}{1}\n'.format("Job name",settings.name))
	log.write('{0:50}{1}\n'.format("Structure input",settings.ircfile))
	log.write('{0:50}{1}\n'.format("Structure input filetype",settings.filetype))
	if settings.analysis:
		log.write('{0:50}{1}\n'.format("\nAnalysis will be performed\n","YES"))
		if len(settings.geo_dist) > 0:
			log.write('    {0:50}\n'.format("The following bond distances will be examiend:\n"))
			for i in range(len(settings.geo_dist)):
				log.write('         {0:15}atom 1: {1}\n         {2:15}atom 2: {3}\n\n'.format("Bond "+str(i+1),settings.geo_dist[i][0],"",settings.geo_dist[i][1]))
		if len(settings.geo_ang) > 0:
			log.write('    {0:50}\n'.format("The following angles will be examiend:\n"))
			for i in range(len(settings.geo_ang)):
				log.write('         {0:15}atom 1: {1}\n         {2:15}atom 2: {3}\n         {2:15}atom 3: {4}\n\n'.format("Angle "+str(i+1),settings.geo_ang[i][0],"",settings.geo_ang[i][1],settings.geo_ang[i][2]))
		if len(settings.geo_dih) > 0:
			log.write('    {0:50}\n'.format("The following dihedral angles will be examiend:\n"))
			for i in range(len(settings.geo_dih)):
				log.write('         {0:15}atom 1: {1}\n         {2:15}atom 2: {3}\n         {2:15}atom 3: {4}\n         {2:15}atom 4: {5}\n\n'.format("Dihedral "+str(i+1),settings.geo_dih[i][0],"",settings.geo_dih[i][1],settings.geo_dih[i][2],settings.geo_dih[i][3]))	
	log.write('{0:50}{1}\n'.format("Keep xyz files","YES" if settings.keepxyz else "NO"))
	log.write('{0:50}{1}\n'.format("Keep input files","YES" if settings.keepinput else "NO"))
	log.write('{0:50}{1}\n'.format("Keep output files","YES" if settings.keepoutput else "NO"))
	log.write('{0:50}{1}\n'.format("Keep this log file","YES" if settings.keeplog else "NO"))
	log.write('{0:50}{1}\n'.format("Reorder structures","YES" if settings.reorder else "NO"))
	log.write('{0:50}{1}\n'.format("Reduce structures","YES" if settings.reduce else "NO"))
	if settings.reduce:
		log.write('{0:50}{1}\n'.format("Threshold for reducing structures",settings.reduce_tresh))
	log.write('{0:50}{1}\n'.format("Prepare input and exit","YES" if settings.prepareonly else "NO"))
	# structure fragment details
	
	
	
	log.close()
	
#-------------------------------------Reorder and reduce----------------------------
	
	
#----------------------------------------------------------------------------------------
	
	
	
	

#-------------------------------------Clean up-----------------------------------------
	
		


if __name__ == "__main__":
    main()