
import numpy as np

# rotates structure, takes strcuture A and rotation R, returns rotated structure B
def rotate(A, rot):
	B=np.dot(A, rot)
	return B

# does partial procrustes analysis, takes structure A and reference R and returns RMSD
def procrustes(A, R):
	rot=find_rotation(A, R)
	#calculates rotated structure B
	B=rotate(A, rot)
	RMSD=calc_RMSD(B,R)
	return RMSD
	
#function that finds centroid from numpy array, returns numpy array with position of centroid
def find_centroid(P):
	C=P.mean(axis=0)
	return C
	
#function that calculates RMSD between two xyz structures, takes structures as (m*3) numpy array
def calc_RMSD(A, B):
	RMSD=0.0
	for i in range(3):
		for x in range(len(A)):
			RMSD += (A[x,i]-B[x,i])**2
	return np.sqrt(RMSD/len(A))

	#takes structure A and reference structure R, calculates optimal rotation (R) using SVD for alignment based on Kabsch algorithm!
def find_rotation(A, R):
	# Computation of the covariance matrix
	C = np.dot(np.transpose(A), R)
	
	u, s, vh = np.linalg.svd(C)
	
	#assure right handed coordinate system)
	if (np.linalg.det(u) * np.linalg.det(vh)) < 0.0:
		s[-1] = -s[-1]
		u[:, -1] = -u[:, -1]
	
	#calculate rotation R
	rot = np.dot(u, vh)
	return rot

def center_xyz(structure, atoms, mode):
	
	#first calculate translation vector V by calculating center and reversing the vector
	if mode == 'c':
		V=-1*find_centroid(structure)
	elif mode == 'm':
		V=-1*find_centerofmass(structure, atoms)
	else:
		sys.exit('something went wrong')
	
	#now add vector to each atom
	structure=structure+V
	return structure

def reorder(structures):
	if len(structures.xyz) == 1: # only one structure
		return structures
	RMSD = []
	for i in range(0, len(structures.xyz)-1):
		RMSD = RMSD + [procrustes(structures.xyz[i],structures.xyz[i+1])]
	if np.amax(RMSD) > 3* np.mean(RMSD):
		cut = RMSD.index(np.amax(RMSD))
		rmsd1=procrustes(structures.xyz[0], structures.xyz[cut+1])
		rmsd2=procrustes(structures.xyz[cut], structures.xyz[cut+1])
		rmsd3=procrustes(structures.xyz[0], structures.xyz[-1])
		rmsd4=procrustes(structures.xyz[cut], structures.xyz[-1])
		
		if rmsd1 < rmsd2 and rmsd1 < rmsd3 and rmsd1 < rmsd4:
			structures.xyz = structures.xyz[cut::-1]+structures.xyz[cut+1:]
			structures.xyz_1 = structures.xyz_1[cut::-1]+structures.xyz_1[cut+1:]
			structures.xyz_2 = structures.xyz_2[cut::-1]+structures.xyz_2[cut+1:]
			structures.title = structures.title[cut::-1]+structures.title[cut+1:]
		elif rmsd2 < rmsd1 and rmsd2 < rmsd3 and rmsd2 < rmsd4:
			pass		
		elif rmsd3 < rmsd1 and rmsd3 < rmsd2 and rmsd3 < rmsd4:
			structures.xyz = structures.xyz[cut::-1]+structures.xyz[:cut:-1]
			structures.xyz_1 = structures.xyz_1[cut::-1]+structures.xyz_1[:cut:-1]
			structures.xyz_2 = structures.xyz_2[cut::-1]+structures.xyz_2[:cut:-1]
			structures.title = structures.title[cut::-1]+structures.title[:cut:-1]
		elif rmsd4 < rmsd1 and rmsd4 < rmsd2 and rmsd4 < rmsd3:
			structures.xyz = structures.xyz[:cut+1]+structures.xyz[:cut:-1]
			structures.xyz_1 = structures.xyz_1[:cut+1]+structures.xyz_1[:cut:-1]
			structures.xyz_2 = structures.xyz_2[:cut+1]+structures.xyz_2[:cut:-1]
			structures.title = structures.title[:cut+1]+structures.title[:cut:-1]
		else:
			pass
	return structures
	
def reduce(structures, settings):
	if len(structures.xyz) == 1: # only one structure
		return structures
	i=0
	count = 0
	while i < len(structures.xyz)-1:
		if procrustes(structures.xyz[i], structures.xyz[i+1]) >= settings.reduce_tresh:
			i = i+1
		else:
			count = count +1
			del structures.xyz[i+1]
			del structures.title[i+1]
			del structures.xyz_1[i+1]
			del structures.xyz_2[i+1]
	log = open(settings.logfile, 'a')
	log.write("Reduced structures by {0} ({1:.0f}%)!\n".format(count, float(count)/(len(structures.xyz)+count)*100))
	log.close()
	return structures