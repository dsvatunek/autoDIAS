#!/usr/bin/env python 

banner="""
*****************************************************************************************
                                   autoDIAS
								Dennis Svatunek
								    UCLA
							d.svatunek@chem.ucla.com
							    @clickchemist
							
						Download the newest version from:
					   https://github.com/dsvatunek/autoDIAS
					   
					 Please cite. DOI: 10.1002/jcc.26023
*****************************************************************************************				   

Python wrapper for several computational chemistry software packages to automatically 
perform a Distortion/Interaction Activation Strain analysis calculations.

*****************************************************************************************
"""

from DIAS_inputparser import *
from DIAS_analysis import *
from DIAS_reduce_reorder import *
import argparse, sys, time, subprocess, os, datetime, copy

__version__= '1.0.0'
__author__= 'Dennis Svatunek'
__email__ = "d.svatunek@chem.ucla.edu"

#Python2 compatibility
try:
	range = xrange
except NameError:
	pass

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
	coordinates = ""
	for i in range(len(A)):
		coordinates = coordinates +("{0:2s} {1:15.12f} {2:15.12f} {3:15.12f}\n".format(atoms[i], A[i, 0], A[i, 1], A[i, 2]))
		#use G09 input string ans substitute the correct things
	g09structure = g09structure.replace("$multiplicity", str(multiplicity))
	g09structure = g09structure.replace("$charge", str(charge))
	g09structure = g09structure.replace("$filename", jobname)
	g09structure = g09structure.replace("$coordinates", coordinates.strip())
	file.write(g09structure)
	file.close()
	return

def molecular_formula(atoms):
	alphabetic_periodic_table = sorted(copy.deepcopy(periodic_table))
	formula = ""
	count = []
	for element in alphabetic_periodic_table:
		count =  count + [atoms.count(element)]
	
	if count[alphabetic_periodic_table.index("C")] > 1:
		formula = formula + 'C' + str(count[alphabetic_periodic_table.index("C")])
		count[alphabetic_periodic_table.index("C")] = 0
	elif count[alphabetic_periodic_table.index("C")] > 0:
		formula = formula + 'C'
		count[alphabetic_periodic_table.index("C")] = 0
	
	if count[alphabetic_periodic_table.index("H")] > 1:
		formula = formula + 'H' + str(count[alphabetic_periodic_table.index("H")])
		count[alphabetic_periodic_table.index("H")] = 0	
	elif count[alphabetic_periodic_table.index("H")] > 0:
		formula = formula + 'H'
		count[alphabetic_periodic_table.index("H")] = 0
		
	for x in range(0, len(alphabetic_periodic_table)):
		if count[x] == 1:
			formula = formula + alphabetic_periodic_table[x]
		elif count[x] > 1:
			formula = formula + alphabetic_periodic_table[x] + str(count[x])
	return formula

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
	settings, structures = parse_in(args.inp_file,args.analysis)
	#analysis only run analysis only and then stop
	if args.analysis:
		endtime=time.time()
		totaltime=str(endtime-starttime)
		seconds=totaltime.split('.')[0]
		milliseconds=float('0.'+totaltime.split('.')[1])*1000	
		print('Processed input and parsed structures in {0} seconds and {1:.0f} ms'.format(seconds, float(milliseconds)))	
		if settings.reorder:
			structures = reorder(structures)
		if settings.reduce:
			structures = reduce(structures, settings)
		analysis_only(structures,settings)
		return
		
	endtime=time.time()
	totaltime=str(endtime-starttime)
	seconds=totaltime.split('.')[0]
	milliseconds=float('0.'+totaltime.split('.')[1])*1000
	
	log = open(settings.logfile, 'a')
	log.write(banner)
	try: #python2.6
		hostname = subprocess.check_output(['hostname']).decode('utf-8')
		username = subprocess.check_output(['whoami']).decode('utf-8')
	except AttributeError:
		hostname= ""
		username= ""
	log.write(str(username).strip()+'@'+str(hostname).strip()+'\n')
	log.write(time.ctime()+'\n')
	log.write('-------------------------------------------------------\n')
	log.write('Input file: '+args.inp_file)
	log.write('\nProcessed input and parsed structures in {0} seconds and {1:.0f} ms\n'.format(seconds, float(milliseconds)))
	log.write('-------------------------------------------------------\n')
	log.write('Determined settings and variables: \n\n')
	log.write('{0:50}{1}\n'.format("Job name",settings.name))
	log.write('{0:50}{1}\n'.format("Structure input",settings.ircfile))
	log.write('{0:50}{1}\n'.format("Structure input filetype",settings.filetype))
	if settings.analysis:
		log.write('{0:50}{1}\n'.format("\nAnalysis will be performed","YES"))
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
	formula = molecular_formula(structures.atoms)
	fragment1_formula = molecular_formula(structures.frag1atoms)
	fragment2_formula = molecular_formula(structures.frag2atoms)
	
	log.write('\nStructure information\n\n')
	log.write('Input structure has {0} atoms with a molecular formula of:{1:>15}\n\n'.format(len(structures.atoms), formula))
	log.write('{0} has {1} atoms with a molecular formula of:{2:>15}\n\n'.format(settings.frag1name, len(structures.frag1atoms), fragment1_formula))
	log.write('{0} has {1} atoms with a molecular formula of:{2:>15}\n\n'.format(settings.frag2name, len(structures.frag2atoms), fragment2_formula))
	log.write('Complex has a charge and multiplicity of: {0} {1}\n'.format(settings.charge, settings.multi))
	log.write('{0} has a charge and multiplicity of: {1} {2}\n'.format(settings.frag1name, settings.frag1charge, settings.frag1multi))
	log.write('{0} has a charge and multiplicity of: {1} {2}\n'.format(settings.frag2name, settings.frag2charge, settings.frag2multi))
	log.write('{0} has an energy of: {1}\n'.format(settings.frag1name, settings.frag1energy))
	log.write('{0} has an energy of: {1}\n'.format(settings.frag2name, settings.frag2energy))
	log.write('\n')
	
	log.write('-------------------------------------------------------\n')		
	log.close()
	
#-------------------------------------Reorder and reduce----------------------------
	if settings.reorder:
		structures = reorder(structures)
	if settings.reduce:
		structures = reduce(structures, settings)
#-------------------------------------Produce input---------------------------------
	starttime=time.time()

	if settings.keepxyz:
		if not os.path.exists(settings.name+'_xyz'):
			os.makedirs(settings.name+'_xyz')
		for x in range(0, len(structures.xyz)):
			save_xyz(structures.xyz[x], structures.atoms, structures.title[x], settings.name+'_xyz/complex_{0:04d}'.format(x+1))
			save_xyz(structures.xyz_1[x], structures.frag1atoms,settings.frag1name+'_{0:04d}'.format(x+1), settings.name+'_xyz/'+settings.frag1name+'_{0:04d}'.format(x+1))
			save_xyz(structures.xyz_2[x], structures.frag2atoms,settings.frag2name+'_{0:04d}'.format(x+1), settings.name+'_xyz/'+settings.frag2name+'_{0:04d}'.format(x+1))
		os.system("cat {0}_xyz/complex*.xyz > {0}_xyz/complete.xyz".format(settings.name))
	if not os.path.exists(settings.name+'_input'):
		os.makedirs(settings.name+'_input')
	if not os.path.exists(settings.name+'_output'):
		os.makedirs(settings.name+'_output')
	
	for x in range (0, len(structures.xyz)):
		save_input(structures.xyz[x], structures.atoms, structures.title[x], settings.name+'_input/complex_{0:04d}.'.format(x+1)+settings.input_file_extension, settings.name+'_output/complex_{0:04d}'.format(x+1), settings.inputlayout, settings.charge, settings.multi)
		save_input(structures.xyz_1[x], structures.frag1atoms, settings.frag1name+'_{0:04d}'.format(x+1), settings.name+'_input/'+settings.frag1name+'_{0:04d}.'.format(x+1)+settings.input_file_extension, settings.name+'_output/'+settings.frag1name+'_{0:04d}'.format(x+1), settings.inputlayout, settings.frag1charge, settings.frag1multi)
		save_input(structures.xyz_2[x], structures.frag2atoms, settings.frag2name+'_{0:04d}'.format(x+1), settings.name+'_input/'+settings.frag2name+'_{0:04d}.'.format(x+1)+settings.input_file_extension, settings.name+'_output/'+settings.frag2name+'_{0:04d}'.format(x+1), settings.inputlayout, settings.frag2charge, settings.frag2multi)	
	
	endtime=time.time()	
	totaltime=str(endtime-starttime)
	seconds=totaltime.split('.')[0]
	milliseconds=float('0.'+totaltime.split('.')[1])*1000
	log = open(settings.logfile, 'a')
	log.write('Produced input files for {0} structures in {1} seconds and {2:.0f} ms\n'.format(len(structures.xyz),seconds, float(milliseconds)))
	log.write('-------------------------------------------------------\n')	
	log.close()	
	
#-------------------------------------Run jobs-----------------------------------------
	if settings.analysis and not settings.prepareonly:
		create_analysis_file(settings)
	log = open(settings.logfile, 'a')
	
	if settings.prepareonly:
		log.write('No calculations requested!\n')
	else:
		for x in range (0, len(structures.xyz)):
			#produce string with correct bash command
			command = settings.submit_setting.replace("$input", settings.name+'_input/complex_{0:04d}.'.format(x+1)+settings.input_file_extension)
			command = command.replace("$output", settings.name+'_output/complex_{0:04d}.'.format(x+1)+settings.output_file_extension)
			log.write(time.ctime()+'{0:<30}'.format('\tStarting {0:>15}'.format('complex_{0:04d}'.format(x+1))))
			log.flush()	
			starttime=time.time()
			os.system(command)
			endtime=time.time()
			totaltime=endtime-starttime
			log.write('\t finished after {0:0>15}\n'.format(str(datetime.timedelta(seconds=totaltime))))	
			log.flush()
			
			command = settings.submit_setting.replace("$input", settings.name+'_input/'+settings.frag1name+'_{0:04d}.'.format(x+1)+settings.input_file_extension)
			command = command.replace("$output", settings.name+'_output/'+settings.frag1name+'_{0:04d}.'.format(x+1)+settings.output_file_extension)
			log.write(time.ctime()+'{0:<30}'.format('\tStarting {0:>15}'.format(settings.frag1name+'_{0:04d}'.format(x+1))))
			log.flush()	
			starttime=time.time()
			os.system(command)
			endtime=time.time()
			totaltime=endtime-starttime
			log.write('\t finished after {0:0>15}\n'.format(str(datetime.timedelta(seconds=totaltime))))	
			log.flush()
			
			command = settings.submit_setting.replace("$input", settings.name+'_input/'+settings.frag2name+'_{0:04d}.'.format(x+1)+settings.input_file_extension)
			command = command.replace("$output", settings.name+'_output/'+settings.frag2name+'_{0:04d}.'.format(x+1)+settings.output_file_extension)
			log.write(time.ctime()+'{0:<30}'.format('\tStarting {0:>15}'.format(settings.frag2name+'_{0:04d}'.format(x+1))))
			log.flush()	
			starttime=time.time()
			os.system(command)
			endtime=time.time()
			totaltime=endtime-starttime
			log.write('\t finished after {0:0>15}\n'.format(str(datetime.timedelta(seconds=totaltime))))	
			log.flush()
			if settings.analysis:
				analysis(settings, structures, x)
				
	totaltime=time.time()-overallstarttime
	log.write('-------------------------------------------------------\n')
	log.write(time.ctime()+' Finished after {0:0>15}\n'.format(str(datetime.timedelta(seconds=totaltime))))
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
	
