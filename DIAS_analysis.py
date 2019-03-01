radian2degree= 57.2958
hartree2kcal = 627.509

import numpy as np
import os, time

def getSCF(file):
	SCF = getmultiSCF(file)

	return SCF

def getmultiSCF(g09out):
	out = open(g09out, "r")
	SCFenergy = 0.00
	out.seek(0) #returns to top of .log file
	out = reversed(list(out)) # reverse list to read from back to find last SCF optimization
	for line in out:
		if "SCF Done" in line: # Gaussian
			SCFenergy = line.split()[4]
			break
		if "Energy=" in line: # Gaussian Force Field
			SCFenergy = line.split()[3]
			break
		elif "FINAL SINGLE POINT ENERGY" in line: # Orca
			SCFenergy = line.split()[4]
			break
		elif "Total energy in the final basis set" in line: # Qchem
			SCFenergy = line.split()[8]
			break
		elif "Total DFT energy" in line: # NWChem, but only for DFT
			SCFenergy = line.split()[4]
			break	
	#check if it is indeed a float
	try:
		SCFenergy = float(SCFenergy)
		pass
	except ValueError:
		SCFenergy = 0.00
	return SCFenergy

def getDistance(A, B):
	distance = np.sqrt( (A[0]-B[0])**2 + (A[1]-B[1])**2 + (A[2]-B[2])**2)
	return distance
	
# get angle between three atoms A, B and C, takes numpy arrays for coordinates
def getAngle(A, B, C):
	angle = 0.00
	#get unit vectors of vector BA and BC
	BA_u = (A - B) / np.linalg.norm(A-B)
	BC_u = (C - B) / np.linalg.norm(C-B)
	#calculate angle
	angle = np.arccos(np.clip(np.dot(BA_u, BC_u), -1, 1))
	angle = angle * radian2degree
	return angle

def getDihedral(p0, p1, p2, p3):
    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

def create_analysis_file(settings):

	analysis_file = open(settings.name+'_DIAS.txt', 'w')
	analysis_file.write('{0:<15}'.format('Step'))
	
	try:
		for x in range(0, len(settings.geo_dist)):
			analysis_file.write('{0:<18}'.format('Distance between'))
	except NameError:
		pass
		
	try:
		for x in range(0, len(settings.geo_ang)):
			analysis_file.write('{0:<18}'.format('Angle between'))
	except NameError:
		pass	
		
	try:
		for x in range(0, len(settings.geo_dih)):
			analysis_file.write('{0:<18}'.format('Dihedral between'))
	except NameError:
		pass	

	
	analysis_file.write('{0:<15}'.format('Total energy'))
	analysis_file.write('{0:<15}'.format('E(Int)'))
	analysis_file.write('{0:<15}'.format('E(Dist)'))	
	analysis_file.write('{0:<15}'.format('E(Dist)'))
	analysis_file.write('{0:<15}'.format('E(Dist)'))
	analysis_file.write('{0:<20}'.format('E(SCF)'))
	analysis_file.write('{0:<20}'.format('E(SCF)'))
	analysis_file.write('{0:<20}'.format('E(SCF)'))
	analysis_file.write('\n')


	analysis_file.write('{0:<15}'.format(''))
	try:
		for x in range(0, len(settings.geo_dist)):
			analysis_file.write('{0:<18}'.format(str(settings.geo_dist[x][0])+'-'+str(settings.geo_dist[x][1])))
	except NameError:
		pass		
	try:
		for x in range(0, len(settings.geo_ang)):
			analysis_file.write('{0:<18}'.format(settings.geo_ang[x][0]+'-'+settings.geo_ang[x][1]+'-'+ settings.geo_ang[x][2]))
	except NameError:
		pass
	try:
		for x in range(0, len(settings.geo_dih)):
			analysis_file.write('{0:<18}'.format(settings.geo_dih[x][0]+'-'+settings.geo_dih[x][1]+'-'+ settings.geo_dih[x][2]+'-'+ settings.geo_dih[x][3]))
	except NameError:
		pass

	analysis_file.write('{0:<15}'.format(''))
	analysis_file.write('{0:<15}'.format(''))
	analysis_file.write('{0:<15}'.format('Total'))	
	analysis_file.write('{0:<15}'.format(settings.frag1name))
	analysis_file.write('{0:<15}'.format(settings.frag2name))
	analysis_file.write('{0:<20}'.format('complex'))
	analysis_file.write('{0:<20}'.format(settings.frag1name))
	analysis_file.write('{0:<20}'.format(settings.frag2name))
	analysis_file.write('\n')

	analysis_file.write('{0:<15}'.format(''))
	try:
		for x in range(0, len(settings.geo_dist)):
			analysis_file.write('{0:<18}'.format('angstrom'))
	except NameError:
		pass			
	try:
		for x in range(0, len(settings.geo_ang)):
			analysis_file.write('{0:<18}'.format('degree'))
	except NameError:
		pass
	try:
		for x in range(0, len(settings.geo_dih)):
			analysis_file.write('{0:<18}'.format('degree'))
	except NameError:
		pass
		
	analysis_file.write('{0:<15}'.format('kcal/mol'))
	analysis_file.write('{0:<15}'.format('kcal/mol'))
	analysis_file.write('{0:<15}'.format('kcal/mol'))	
	analysis_file.write('{0:<15}'.format('kcal/mol'))
	analysis_file.write('{0:<15}'.format('kcal/mol'))
	analysis_file.write('{0:<20}'.format('hartree'))
	analysis_file.write('{0:<20}'.format('hartree'))
	analysis_file.write('{0:<20}'.format('hartree'))
	analysis_file.write('\n')
	try:
		for x in range(0, len(settings.geo_dist)):
			analysis_file.write('{0:-<18}'.format(''))
	except NameError:
		pass
	try:
		for x in range(0, len(settings.geo_ang)):
			analysis_file.write('{0:-<18}'.format(''))
	except NameError:
		pass
	try:
		for x in range(0, len(settings.geo_dih)):
			analysis_file.write('{0:-<18}'.format(''))
	except NameError:
		pass
	analysis_file.write('{0:-<150}'.format(''))
	analysis_file.write('\n')
	analysis_file.close()
	return

def analysis(settings, structures, x):
	analysis_file = open(settings.name+'_DIAS.txt', 'a')
	irc_SCF = getSCF(settings.name+'_output/complex_{0:04d}.'.format(x+1)+settings.output_file_extension)
	fragment1_SCF = getSCF(settings.name+'_output/'+settings.frag1name+'_{0:04d}.'.format(x+1)+settings.output_file_extension)
	fragment2_SCF = getSCF(settings.name+'_output/'+settings.frag2name+'_{0:04d}.'.format(x+1)+settings.output_file_extension)

	# calculate distortion energiesdef analysis_only(structures,settings):
	fragment1_dist = fragment1_SCF - settings.frag1energy	
	fragment2_dist = fragment2_SCF - settings.frag2energy
	# calculate interaction energy
	interaction_energy = irc_SCF - (fragment1_SCF + fragment2_SCF)

	# extract geometric values
	distance_values = []
	for y in range(0, len(settings.geo_dist)):
		distance_values = distance_values + [getDistance(structures.xyz[x][int(settings.geo_dist[y][0])-1],structures.xyz[x][int(settings.geo_dist[y][1])-1])]

	angle_values = []
	for y in range(0, len(settings.geo_ang)):
		angle_values = angle_values + [getAngle(structures.xyz[x][int(settings.geo_ang[y][0])-1],structures.xyz[x][int(settings.geo_ang[y][1])-1],structures.xyz[x][int(settings.geo_ang[y][2])-1])]	# -1 since array begins with 0
	dihedral_values = []
	for y in range(0, len(settings.geo_dih)):
		dihedral_values = dihedral_values + [getDihedral(structures.xyz[x][int(settings.geo_dih[y][0])-1],structures.xyz[x][int(settings.geo_dih[y][1])-1],structures.xyz[x][int(settings.geo_dih[y][2])-1],structures.xyz[x][int(settings.geo_dih[y][3])-1])]	# -1 since array begins with 0

	# print to file
	analysis_file.write('{0:<15}'.format('{0:04d}'.format(x+1)))

	for element in distance_values:
		analysis_file.write('{0:<18}'.format('{0:.5f}'.format(float(element))))
	for element in angle_values:
		analysis_file.write('{0:<18}'.format('{0:.3f}'.format(float(element))))
	for element in dihedral_values:
		analysis_file.write('{0:<18}'.format('{0:.3f}'.format(float(element))))

	analysis_file.write('{0:<15}'.format('{0:.5f}'.format(float((irc_SCF-settings.frag1energy-settings.frag2energy)*hartree2kcal))))
	analysis_file.write('{0:<15}'.format('{0:.5f}'.format(float(interaction_energy*hartree2kcal))))
	analysis_file.write('{0:<15}'.format('{0:.5f}'.format(float((fragment1_dist+fragment2_dist)*hartree2kcal))))
	analysis_file.write('{0:<15}'.format('{0:.5f}'.format(float(fragment1_dist*hartree2kcal))))
	analysis_file.write('{0:<15}'.format('{0:.5f}'.format(float(fragment2_dist*hartree2kcal))))
	analysis_file.write('{0:<20}'.format('{0:.9f}'.format(float(irc_SCF))))
	analysis_file.write('{0:<20}'.format('{0:.9f}'.format(float(fragment1_SCF))))
	analysis_file.write('{0:<20}'.format('{0:.9f}'.format(float(fragment2_SCF))))
	analysis_file.write('\n')		
	analysis_file.close()
	
	
def analysis_only(structures, settings):
	starttime= time.time()
	create_analysis_file(settings)
	for i in range(0, len(structures.xyz)):
		analysis(settings, structures, i)
	os.remove(settings.logfile)
	endtime=time.time()
	totaltime=str(endtime-starttime)
	seconds=totaltime.split('.')[0]
	milliseconds=float('0.'+totaltime.split('.')[1])*1000	
	print('Analysis of {0} structures done in {1} seconds and {2:.0f} ms'.format(len(structures.xyz),seconds, float(milliseconds)))	
	return