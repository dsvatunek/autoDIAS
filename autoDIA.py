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
from DIA_analysis import *
import argparse, sys, time, subprocess, os, datetime

__version__= ' '

#Python2 compatibility
try:
	range = xrange
except NameError:
	pass

radian2degree= 57.2958
hartree2kcal = 627.509

def save_xyz(A, atoms, title, filename):
	file = open(filename + '.xyz', 'w')
	file.write(str(len(atoms))+ "\n"+ title + "\n")
	for i in range(len(A)):
		file.write("{0:2s} {1:15.12f} {2:15.12f} {3:15.12f}\n".format(atoms[i], A[i, 0], A[i, 1], A[i, 2]))
	file.close()
	return

def save_input(A, atoms, title, filename, jobname, g09structure, charge, multiplicity):
	file = open(filename, 'w')
	#make a string with coordinates
	coordinates = str(charge) + ' ' + str(multiplicity) + '\n'
	for i in range(len(A)):
		coordinates = coordinates +("{0:2s} {1:15.12f} {2:15.12f} {3:15.12f}\n".format(atoms[i], A[i, 0], A[i, 1], A[i, 2]))
		#use G09 input string ans substitute the correct things
	g09structure = g09structure.replace("$filename", jobname)
	g09structure = g09structure.replace("$coordinates", coordinates.strip())
	file.write(g09structure)
	file.close()
	return


def main():

	parser =  argparse.ArgumentParser(usage='%(prog)s inp_file')
	parser.add_argument('inp_file', metavar='Input file', type=str, help='Name of the input file')
	parser.add_argument("-a", "--analysis",  action='store_true', help='Analysis only', default=False)
	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit(1)	
	args = parser.parse_args()
	
#----------------------------------------------------------------------------------------
	overallstarttime=time.time()	
#-------------------------------------Parse Settings and Inputs----------------------
	starttime=time.time() # start timing for input parsing
	settings, structures = parse_in(args.inp_file)
	#analysis only run analysis only and then stop
	if args.analysis:
		print("analysis only")
		analysis_only(structures,settings)
		return
		
		
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
			log.write('    {0:50}\n'.format("The following bond distances will be examined:\n"))
			for i in range(len(settings.geo_dist)):
				log.write('         {0:15}atom 1: {1}\n         {2:15}atom 2: {3}\n\n'.format("Bond "+str(i+1),settings.geo_dist[i][0],"",settings.geo_dist[i][1]))
		if len(settings.geo_ang) > 0:
			log.write('    {0:50}\n'.format("The following angles will be examined:\n"))
			for i in range(len(settings.geo_ang)):
				log.write('         {0:15}atom 1: {1}\n         {2:15}atom 2: {3}\n         {2:15}atom 3: {4}\n\n'.format("Angle "+str(i+1),settings.geo_ang[i][0],"",settings.geo_ang[i][1],settings.geo_ang[i][2]))
		if len(settings.geo_dih) > 0:
			log.write('    {0:50}\n'.format("The following dihedral angles will be examined:\n"))
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
	
	
	log.write('-------------------------------------------------------\n')		
	log.close()
	
#-------------------------------------Reorder and reduce----------------------------
	
	
#-------------------------------------Produce input---------------------------------
	starttime=time.time()
	
	#get fragment xyz
	structures.xyz_1 = []
	settings.frag1atoms[:] = [x -1 for x in settings.frag1atoms]
	for x in range(0, len(structures.xyz)):
		structures.xyz_1 = structures.xyz_1 + [structures.xyz[x][settings.frag1atoms]]
	structures.xyz_2 = []
	settings.frag2atoms[:] = [x -1 for x in settings.frag2atoms]
	for x in range(0, len(structures.xyz)):
		structures.xyz_2 = structures.xyz_2 + [structures.xyz[x][settings.frag2atoms]]
	structures.frag1atoms = []
	for element in settings.frag1atoms:
		structures.frag1atoms = structures.frag1atoms + [structures.atoms[element]]
	structures.frag2atoms = []
	for element in settings.frag2atoms:
		structures.frag2atoms = structures.frag2atoms + [structures.atoms[element]]

	if settings.keepxyz:
		if not os.path.exists(settings.name+'_xyz'):
			os.makedirs(settings.name+'_xyz')
		for x in range(0, len(structures.xyz)):
			save_xyz(structures.xyz[x], structures.atoms, structures.title[x], settings.name+'_xyz/irc_{0:04d}'.format(x+1))
			save_xyz(structures.xyz_1[x], structures.frag1atoms,settings.frag1name+'_{0:04d}'.format(x+1), settings.name+'_xyz/'+settings.frag1name+'_{0:04d}'.format(x+1))
			save_xyz(structures.xyz_2[x], structures.frag2atoms,settings.frag2name+'_{0:04d}'.format(x+1), settings.name+'_xyz/'+settings.frag2name+'_{0:04d}'.format(x+1))		
	if not os.path.exists(settings.name+'_input'):
		os.makedirs(settings.name+'_input')
	if not os.path.exists(settings.name+'_output'):
		os.makedirs(settings.name+'_output')
	
	for x in range (0, len(structures.xyz)):
		save_input(structures.xyz[x], structures.atoms, structures.title[x], settings.name+'_input/irc_{0:04d}.'.format(x+1)+settings.input_file_extension, settings.name+'_output/irc_{0:04d}'.format(x+1), settings.inputlayout, settings.charge, settings.multi)
		save_input(structures.xyz_1[x], structures.frag1atoms, settings.frag1name+'_{0:04d}'.format(x+1), settings.name+'_input/'+settings.frag1name+'_{0:04d}.'.format(x+1)+settings.input_file_extension, settings.name+'_output/'+settings.frag1name+'_{0:04d}'.format(x+1), settings.inputlayout, settings.frag1charge, settings.frag1multi)
		save_input(structures.xyz_2[x], structures.frag2atoms, settings.frag2name+'_{0:04d}'.format(x+1), settings.name+'_input/'+settings.frag2name+'_{0:04d}.'.format(x+1)+settings.input_file_extension, settings.name+'_output/'+settings.frag2name+'_{0:04d}'.format(x+1), settings.inputlayout, settings.frag2charge, settings.frag2multi)	
	
	endtime=time.time()	
	totaltime=str(endtime-starttime)
	seconds=totaltime.split('.')[0]
	milliseconds=float('0.'+totaltime.split('.')[1])*1000
	log = open(settings.logfile, 'a')
	log.write('\nProduced input files in {} seconds and {:.0f} ms\n'.format(seconds, float(milliseconds)))
	log.write('-------------------------------------------------------\n')	
	log.close()	
	
#-------------------------------------Run jobs-----------------------------------------
	create_analysis_file(settings)
	log = open(settings.logfile, 'a')
	
	if settings.prepareonly:
		log.write('No calculations requested!')
	else:
		for x in range (0, len(structures.xyz)):
			#produce string with correct bash command
			command = settings.submit_setting.replace("$input", settings.name+'_input/irc_{0:04d}.'.format(x+1)+settings.input_file_extension)
			command = command.replace("$output", settings.name+'_output/irc_{0:04d}.'.format(x+1)+settings.output_file_extension)
			log.write(time.ctime()+'{:<30}'.format('\tStarting {:>15}'.format('irc_{0:04d}'.format(x+1))))
			log.flush()	
			starttime=time.time()
			os.system(command)
			endtime=time.time()
			totaltime=endtime-starttime
			log.write('\t finished after {:0>15}\n'.format(str(datetime.timedelta(seconds=totaltime))))	
			log.flush()
			
			command = settings.submit_setting.replace("$input", settings.name+'_input/'+settings.frag1name+'_{0:04d}.'.format(x+1)+settings.input_file_extension)
			command = command.replace("$output", settings.name+'_output/'+settings.frag1name+'_{0:04d}.'.format(x+1)+settings.output_file_extension)
			log.write(time.ctime()+'{:<30}'.format('\tStarting {:>15}'.format(settings.frag1name+'_{0:04d}'.format(x+1))))
			log.flush()	
			starttime=time.time()
			os.system(command)
			endtime=time.time()
			totaltime=endtime-starttime
			log.write('\t finished after {:0>15}\n'.format(str(datetime.timedelta(seconds=totaltime))))	
			log.flush()
			
			command = settings.submit_setting.replace("$input", settings.name+'_input/'+settings.frag2name+'_{0:04d}.'.format(x+1)+settings.input_file_extension)
			command = command.replace("$output", settings.name+'_output/'+settings.frag2name+'_{0:04d}.'.format(x+1)+settings.output_file_extension)
			log.write(time.ctime()+'{:<30}'.format('\tStarting {:>15}'.format(settings.frag2name+'_{0:04d}'.format(x+1))))
			log.flush()	
			starttime=time.time()
			os.system(command)
			endtime=time.time()
			totaltime=endtime-starttime
			log.write('\t finished after {:0>15}\n'.format(str(datetime.timedelta(seconds=totaltime))))	
			log.flush()
			if analysis:
				analysis(settings, structures, x)
				
		totaltime=time.time()-overallstarttime
	log.write('\n\n'+time.ctime()+' Finished after {:0>15}\n'.format(str(datetime.timedelta(seconds=totaltime))))
	log.write('-------------------------------------------------------\n')
	log.write('-------------------------------------------------------\n')
	log.close()	
#-------------------------------------Clean up-----------------------------------------
	if settings.keepinput == False:
		shutil.rmtree("input/", ignore_errors=True)
	if settings.keepoutput == False:
		shutil.rmtree("output/", ignore_errors=True)	
	if settings.keeplog == False:
		os.remove(settings.logfile)
		
		
if __name__ == "__main__":


	main()
	
