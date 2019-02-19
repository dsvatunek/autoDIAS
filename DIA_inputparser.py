import numpy as np
import sys
import re

def line_prepender(filename, line):
	with open(filename, 'r+') as f:
		content = f.read()
		f.seek(0, 0)
		f.write(line.rstrip('\r\n') + '\n' + content)

class structures:
	pass

periodic_table = ["","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
    "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl",
    "Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"]
	
#Covalent radii taken from DOI: 10.1039/b801115j
#Everything beyond Cm was set to 1.80
covalent_radii = [0.00,0.32,0.28,1.28,0.96,0.84,0.76,0.71,0.66,0.57,0.58,1.66,1.41,1.21,1.11,1.07,1.05,1.02,1.06,2.03,1.76,1.70,1.60,1.53,1.39,1.61,1.52,1.50,1.24,1.32,1.22,1.22,1.20,1.19,1.20,1.20,1.16,2.20,1.95,1.90,1.75,
    1.64,1.54,1.47,1.46,1.42,1.39,1.45,1.44,1.42,1.39,1.39,1.38,1.39,1.40,2.44,2.15,2.07,2.04,2.03,2.01,1.99,1.98,1.98,1.96,1.94,1.92,1.92,1.89,1.90,1.87,1.87,1.75,1.70,1.62,1.51,1.44,1.41,1.36,1.36,1.32,1.45,
    1.46,1.48,1.40,1.50,1.50,2.60,2.21,2.15,206,2.00,1.96,1.90,1.87,180,169,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80]

#takes xyz coordinates as numpy array and returns distance
def getDistance(atom1, atom2):
	distance = 0.00
	distance = np.sqrt((atom1[0]-atom2[0])**2+(atom1[1]-atom2[1])**2+(atom1[2]-atom2[2])**2)
	return distance

# checks if "s" is an int, returns true or false
def isInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False	

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
	
def isBond(distance, atomtype1, atomtype2, covalent_radii):
	bond = True
	'''find bond threshold based on covalent radii
	taken from 10.1039/b801115j
	For everything beyond Cm 1.80 A was used.
	get radii from saved numpy array
	use atomtypes as numbers
	use factor to be sure to catch bonds. Bigger factor means more bonds are detected
	'''
	atomtype1 = int(periodic_table.index(atomtype1))
	atomtype2 = int(periodic_table.index(atomtype2))
	bond_threshold = 1.25 * (covalent_radii[atomtype1]+covalent_radii[atomtype2])
	if distance > bond_threshold:
		bond = False
	return bond
	
def getFragments(structure, atoms):
	#get connectivity for each irc structure
	adjacency_matrix = []
	# iterate over every atom 
	for i in range(0, len(structure)):
		distances = []
		#get the distance to every other atom and itself
		for j in range(0, len(structure)):
			distances = distances + [getDistance(structure[i], structure[j])]
		adjacency_matrix = adjacency_matrix + [distances]

	#convert distance to 1 for bond, 0 for no bond
	for i in range(0, len(adjacency_matrix)):
		for j in range(0, len(adjacency_matrix[i])):
			if isBond(adjacency_matrix[i][j], atoms[i], atoms[j], covalent_radii):
				adjacency_matrix[i][j] = 1
			else:
				adjacency_matrix[i][j] = 0
	# now make a list of fragments 
	#first convert each line to a list of the atoms in the structure which show a bond
	for i in range(0, len(adjacency_matrix)):
		adjacency_matrix[i] = [j+1 for j,x in enumerate(adjacency_matrix[i]) if x == 1]
	#make empty fragments list
	fragments = [adjacency_matrix[0]]
	#now iterate over all the elements in adjacency_matrix
	for i in range (1, len(adjacency_matrix)):
		#now iterate over all the fragments in fragments and check if elements are the same
		for j in range(0, len(fragments)):
			#check if elements are already present in this list...
			if len(set(adjacency_matrix[i]) & set(fragments[j])) > 0:
				#....if they are add the lists together and eliminate duplicates...
				fragments[j] = list(set(fragments[j]+adjacency_matrix[i]))
				break
		else:
			#....add this list as new list in fragments
			fragments = fragments + [adjacency_matrix[i]]
	#now we got all of our fragments but some might belong to the same structure, so we need to test if we need to combine them, then combine them, 
	#then check again and so on until nothing else changes. iterative process. Take length of fragments for that since it converged if amount of fragments doesn't changes
	n_frag_old=len(fragments)
	#set j to some value that cannot be the length of fragments
	n_frag_new=-5
	while n_frag_old != n_frag_new:
		#set i to current value of len(fragments)
		n_frag_old=len(fragments)
		new_fragments = [sorted(fragments[0])]
		#combine fragments like before
		for i in range(0, len(fragments)):
			for j in range(0, len(new_fragments)):
				if len(set(fragments[i]) & set(new_fragments[j])) > 0:
					new_fragments[j] = sorted(list(set(new_fragments[j]+fragments[i])))
					break
			else:
				new_fragments = new_fragments + [fragments[i]]
		fragments=new_fragments
		n_frag_new=len(fragments)
	return fragments	

def getConnectivity(structure):
	connectivity_matrix = np.zeros((len(structure),len(structure)))
	for i in range(0,len(connectivity_matrix )):
		for j in range(0,len(connectivity_matrix)):
			distance=getDistance(structure[i],structure[j])
			if isBond(distance,structures.atoms[i], structures.atoms[j],covalent_radii):
				connectivity_matrix[i][j] = 1
			else:
				connectivity_matrix[i][j] = 0
	return connectivity_matrix


# finds up to two bonds that get formed
def getBonds(structures):
	count = np.zeros((structures.xyz[0].shape[0],structures.xyz[0].shape[0]))
	#construct connectivity matrix step 0
	connectivity_old = getConnectivity(structures.xyz[0])	
	#now iterate over all the matrices and count changes
	for element in structures.xyz:
		connectivity = getConnectivity(element)
		for i in range(0, len(connectivity)):
			for j in range(0, len(connectivity[i])):
				if connectivity[i][j] != connectivity_old[i][j]:
					count[i][j] = count[i][j] +1
		connectivity_old= connectivity
	#clear of lower triangle of the matrix to eliminate duplicates
	for i in range(0,len(count)):
		for j in range(0, len(count)):
			if i >= j:
				count[i][j] = 0
	
	count[count == 0] = 1000	
	#now find conenctivity for lowest number in count matrix, equals the bond that changed status the least, might not be the optimal solution
	freq = np.amin(count)
	geodist = [[int(np.where(count == freq)[0][0])+1,int(np.where(count == freq)[1][0])+1]]
	count[np.where(count == freq)[0][0],np.where(count == freq)[1][0]] = 1001
	
	freq = np.amin(count)
	geodist = geodist + [[int(np.where(count == freq)[0][0])+1,int(np.where(count == freq)[1][0])+1]]

	return geodist

	
def duplicates(int_list):
	int_list = sorted(int_list)
	duplicates = False
	for i in range (0, len(int_list)-1):
		if int_list[i] == int_list[i+1]:
			duplicates = True
		else:
		   pass
	return duplicates
	
def auto_frag(structures,settings):
	list_of_fragments = []
	for i in range(0, len(structures.xyz)):
		list_of_fragments =  list_of_fragments + [getFragments(structures.xyz[i], structures.atoms)]
	n_fragments = []
	for i in range(0,len(list_of_fragments)):
		n_fragments = n_fragments + [len(list_of_fragments[i])]
	#this part determines if there is at least one structure with 2 fragments and at least one with more than 2
	twofrag = False
	morefrag = False
	for element in n_fragments:
		if element == 2:
			twofrag = True
			break
		if element > 2:
			morefrag = True
	if twofrag:
		indices = [i for i, x in enumerate(n_fragments) if x == 2]
		counts = [list_of_fragments.count(list_of_fragments[x]) for x in indices]
		fragments = list_of_fragments[indices[counts.index(max(counts))]]
		settings.frag1atoms = np.sort(fragments[0])
		settings.frag2atoms = np.sort(fragments[1])
		print(fragments)
	elif morefrag:
		n_fragments = [x for x in n_fragments if x > 2]
		indices = [i for i, x in enumerate(n_fragments) if x == min(n_fragments)]
		counts = [list_of_fragments.count(list_of_fragments[x]) for x in indices]
		fragments = list_of_fragments[indices[counts.index(max(counts))]]
		settings.frag1atoms = np.sort(fragments[0])
		# now combine all other fragments into one
		settings.frag2atoms = np.sort([j for i in fragments[1:] for j in i])
		line_prepender(settings.logfile, "WARNING: more than two fragments detected! Check fragments!\n\n")
	else:
		line_prepender(settings.logfile, "CRITICAL ERROR: couldn't determine fragments automatically!\n\n")
		sys.exit("CRITICAL ERROR: couldn't determine fragments automatically!")
	return settings
 
def get_single(file):
	structures.atoms = []
	structures.xyz = []
	structures.title = ['single structure']
	with open(file) as input_file:
		file_contents = input_file.read() 
	if "Optimization completed." in file_contents:
		regex = r'\ ---------------------------------------------------------------------\n((?:(?!\ ---------------------------------------------------------------------\n).)*?)\ ---------------------------------------------------------------------\n(?:(?!\ ---------------------------------------------------------------------\n).)*?Stationary\ point\ found'
		for match in re.finditer(regex, file_contents, re.IGNORECASE|re.DOTALL):
			raw_coordinates = match.group(1).strip().split("\n")
			break
	else:
		regex = r'\ ---------------------------------------------------------------------\n((?:(?!\ ---------------------------------------------------------------------\n).)*?)\ ---------------------------------------------------------------------\n(?:(?!\ ---------------------------------------------------------------------\n).)*?Rotational\ constants'
		for match in re.finditer(regex, file_contents, re.IGNORECASE|re.DOTALL):	
			raw_coordinates = match.group(1).strip().split("\n")
		
	structures.xyz.append(np.array([x.split()[3:] for x in raw_coordinates]))
	structures.atoms = [x.split()[1] for x in raw_coordinates]
	structures.atoms = [periodic_table[int(x)] for x in structures.atoms]
	for i in range(len(structures.xyz)):
		structures.xyz[i] = structures.xyz[i].astype(float)
	return structures
	
def get_scan(file):
	with open(file) as input_file:
		file_contents = input_file.read() 
	structures.xyz = []
	structures.atoms = []

	# find everything except the TS structure
	regex = r'---------------------------------------------------------------------((?:(?!---------------------------------------------------------------------).)*?)---------------------------------------------------------------------(?:(?!---------------------------------------------------------------------).)*?Stationary\ point\ found'
	for match in re.finditer(regex, file_contents, re.IGNORECASE|re.DOTALL):
		raw_coordinates = match.group(1).strip().split("\n")
		structures.xyz.append(np.array([x.split()[3:] for x in raw_coordinates]))
	structures.atoms = [x.split()[1] for x in raw_coordinates]
	structures.atoms = [periodic_table[int(x)] for x in structures.atoms]		
		
	for i in range(len(structures.xyz)):
		structures.xyz[i] = structures.xyz[i].astype(float)	
	structures.title = [str(x) for x in range(len(structures.xyz))] 
	return structures

def get_irc(file):
	with open(file) as input_file:
		file_contents = input_file.read() 
	structures.xyz = []
	structures.atoms = []
	# find TS structure, it's the first xyz structure available
	regex = r'Z\n ---------------------------------------------------------------------\n(.*?)\n ---------------------------------------------------------------------\n'
	for match in re.finditer(regex, file_contents, re.IGNORECASE|re.DOTALL):
		raw_coordinates = match.group(1).strip().split("\n")
		break
	structures.xyz.append(np.array([x.split()[3:] for x in raw_coordinates]))
	structures.atoms = [x.split()[1] for x in raw_coordinates]
	structures.atoms = [periodic_table[int(x)] for x in structures.atoms]
	# find everything except the TS structure
	regex = r'---------------------------------------------------------------------((?:(?!---------------------------------------------------------------------).)*?)---------------------------------------------------------------------(?:(?!---------------------------------------------------------------------).)*?Delta-x\ Convergence\ Met'
	for match in re.finditer(regex, file_contents, re.IGNORECASE|re.DOTALL):
		raw_coordinates = match.group(1).strip().split("\n")
		structures.xyz.append(np.array([x.split()[3:] for x in raw_coordinates]))
		
		
	for i in range(len(structures.xyz)):
		structures.xyz[i] = structures.xyz[i].astype(float)	
	structures.title = [str(x) for x in range(len(structures.xyz))]
	return structures
 
def check_fragments(settings,structures):
	for x in settings.frag1atoms:
		if int(x) > len(structures.atoms):
			line_prepender(settings.logfile, "CRITICAL ERROR: an atom number in fragment 1 is higher than the highest atom number!\n\n")
			sys.exit("2")
	for x in settings.frag2atoms:
		if int(x) > len(structures.atoms):
			line_prepender(settings.logfile, "CRITICAL ERROR: an atom number in fragment 2 is higher than the highest atom number!\n\n")
			sys.exit("3")
	if len(settings.frag1atoms) + len(settings.frag2atoms) != len(structures.atoms):
		line_prepender(settings.logfile, "CRITICAL ERROR: number of specified fragment atoms does not equal total number of atoms!\n\n")
		sys.exit("1")
	if duplicates(settings.frag1atoms):
		line_prepender(settings.logfile, "CRITICAL ERROR: at least one atom was defined more than once in fragment1!\n\n")
		sys.exit("4")
	if duplicates(settings.frag2atoms):
		line_prepender(settings.logfile, "CRITICAL ERROR: at least one atom was defined more than once in fragment2!\n\n")
		sys.exit("5")	
	elif len(set(settings.frag1atoms) & set(settings.frag2atoms)) > 0: 
		line_prepender(settings.logfile, "CRITICAL ERROR: at least one atom was defined in both fragments!\n\n")
		sys.exit("6")	
	return

def check_geo(structures, settings):
	natoms = len(structures.xyz[0])
	i = 0
	while i < len(settings.geo_dist): #check if every bond has two elements
		if len(settings.geo_dist[i]) != 2:
			del settings.geo_dist[i]
			line_prepender(settings.logfile, "WARNING: defined bond did not have two elements! - This input is ignored\n\n")
		else:
			i = i+1
	i = 0
	while i < len(settings.geo_dist): #check if it's inly integers
		if isInt(settings.geo_dist[i][0]) and isInt(settings.geo_dist[i][1]):
			i = i+1
		else:
			del settings.geo_dist[i]
			line_prepender(settings.logfile, "WARNING: defined bond did have non-integer input! - This input is ignored\n\n")	
	if len(settings.geo_dist) == 0:
		settings.geo_dist.append(getBonds(structures))
		line_prepender(settings.logfile, "WARNING: no correct definitions of bonds found! - automatic detection of forming bonds\n\n")	
	i = 0
	while i < len(settings.geo_dist): #check if it's out of range
		if int(settings.geo_dist[i][0]) > natoms or int(settings.geo_dist[i][1]) > natoms or int(settings.geo_dist[i][0]) < 1 or int(settings.geo_dist[i][1]) < 1 :
			del settings.geo_dist[i]
			line_prepender(settings.logfile, "WARNING: bond definition: atom number out of range! - This input is ignored\n\n")
		else:
			i = i+1		
	# angles	
	i=0
	while i < len(settings.geo_ang): 
		if len(settings.geo_ang[i]) != 3:
			del settings.geo_ang[i]
			line_prepender(settings.logfile, "WARNING: defined angle did not have three elements! - This input is ignored\n\n")
		else:
			i = i+1
	i = 0
	while i < len(settings.geo_ang): #check if it's inly integers
		if isInt(settings.geo_ang[i][0]) and isInt(settings.geo_ang[i][1]) and isInt(settings.geo_ang[i][2]):
			i = i+1
		else:
			del settings.geo_ang[i]
			line_prepender(settings.logfile, "WARNING: defined angle did have non-integer input! - This input is ignored\n\n")
	i=0
	while i < len(settings.geo_ang): #check if it's out of range
		if int(settings.geo_ang[i][0]) > natoms or int(settings.geo_ang[i][1]) > natoms or int(settings.geo_ang[i][2]) > natoms or int(settings.geo_ang[i][0]) < 1 or int(settings.geo_ang[i][1]) < 1 or int(settings.geo_ang[i][2]) < 1:
			del settings.geo_ang[i]
			line_prepender(settings.logfile, "WARNING: angle definition: atom number out of range! - This input is ignored\n\n")
		else:
			i = i+1					
	# dihedral angles			
	i=0
	while i < len(settings.geo_dih): 
		if len(settings.geo_dih[i]) != 4:
			del settings.geo_dih[i]
			line_prepender(settings.logfile, "WARNING: defined dihedral angle did not have three elements! - This input is ignored\n\n")
		else:
			i = i+1
	i = 0
	while i < len(settings.geo_dih): #check if it's inly integers
		if isInt(settings.geo_dih[i][0]) and isInt(settings.geo_dih[i][1]) and isInt(settings.geo_dih[i][2]) and isInt(settings.geo_dih[i][3]):
			i = i+1
		else:
			del settings.geo_dih[i]
			line_prepender(settings.logfile, "WARNING: defined dihedral angle did have non-integer input! - This input is ignored\n\n")	
	i=0
	while i < len(settings.geo_dih): #check if it's out of range
		if int(settings.geo_dih[i][0]) > natoms or int(settings.geo_dih[i][1]) > natoms or int(settings.geo_dih[i][2]) > natoms or int(settings.geo_dih[i][3]) > natoms or int(settings.geo_dih[i][0]) < 1 or int(settings.geo_dih[i][1]) < 1 or int(settings.geo_dih[i][2]) < 1 or int(settings.geo_dih[i][3]) < 1:
			del settings.geo_dih[i]
			line_prepender(settings.logfile, "WARNING: dihedral angle definition: atom number out of range! - This input is ignored\n\n")
		else:
			i = i+1	
	return settings

def structures_from_G(file, settings):
	file_object = open(file, 'r')
	input = (line for line in file_object)
	type = "single"
	for line in input:
		if "IRC-IRC-IRC" in line:
			type = "IRC"
			break
		elif "Scan" in line:
			type = "scan"
			break
	# Single
	if type == "single":
		structures = get_single(file)
	elif type == "scan":
		structures = get_scan(file)
	else:
		structures = get_irc(file)
	# Get Charge and multiplicity
	file_object.seek(0)
	input_file = (line for line in file_object) # make generator	
	for line in input_file:
		if "Charge =" in line:
			settings.charge = line.split()[2]
			settings.multi = line.split()[5]
			break
	return structures, settings


def structures_from_xyz(file):
	structures.atoms = []
	structures.xyz = []
	structures.title = []
	file_object = open(file, 'r')
	input = (line for line in file_object) #make generator
	#search for number of atoms
	for line in input:
		if isInt(line.strip()):
			n_atoms=int(line)
			break	
	else: #exits if no line with number of atoms was found
		sys.exit('Error:\t\tNo xyz coordinates found in file: ' + file)
	#skip one line
	structures.title.append(next(input).strip())
	# now there should be n_atoms lines of coordinates
	for i in range(n_atoms):
		l=next(input).split()
		if l[0] in periodic_table:
			structures.atoms.append(l[0]) #get atom symbol and append to atom list
		else:
			sys.exit('Error:\t\tsomething is wrong with the first structure in file: '+file)
		coords=[float(x) for x in l[1:]] #convert line to list of floats
		coords=np.array([coords]) #create array with coords
		try: #try append, doesn't work if XYZ doesn't exist yet
			XYZ=np.concatenate((XYZ,coords), axis=0)
		except NameError:
			XYZ=coords
	structures.xyz.append(XYZ) #append first structure to structures list
	del XYZ #get rid of that for the next structure
	#now search for more structures
	for line in input:
		#start extracting if atom number line is found
		try:
			if int(line.strip()) == n_atoms:
				#read one line to skip title
				structures.title.append(next(input).strip())
				# now there should be n_atoms lines of coordinates
				for i in range(n_atoms):
					l=next(input).split()
					coords=[float(x) for x in l[1:]]
					coords=np.array([coords])
					try: #try append, doesn't work if XYZ doesn't exist yet
						XYZ=np.concatenate((XYZ,coords), axis=0)
					except NameError:
						XYZ=coords
				structures.xyz.append(XYZ)
				del XYZ
		except ValueError:
			pass			
	return structures

def parse_in(input_filename, analysisonly):
	class settings:
		pass
	settings.void = []
	input_object = open(input_filename, 'r')
#------------structures filename
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	for line in input_file:
		if len(line.split()) > 1: # This is needed otherwise it produces out of range error when contents of lines too short
			if line.split()[0].upper()+line.split()[1].upper() == "INPUTFILE":
				settings.ircfile = line.split()[-1]
				break
	else: #information on input file type is missing
		print('Error:\t\tNo information on the structure input file could be found')
		sys.exit("7")	
#------------GET JOBNAME
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator
	settings.name = "DIA"
	for line in input_file:
		if len(line.split()) > 1:
			if line.split()[0].upper()+line.split()[1].upper() == "JOBNAME":
				try:
					settings.name = line.split()[-1]
					break
				except IndexError:
					break
#------------Set Logfilename
	if analysisonly:
		settings.logfile= settings.name + "_analysisonly_log.txt"
		f= open(settings.logfile, 'w')
		f.close()
	
	else:	
		settings.logfile= settings.name + "_log.txt"
		f= open(settings.logfile, 'w')
		f.close()
#------------Parse Structures
	#find out what kind of file we are dealing with. See if it's a gaussian file
	try:
		structure_object = open(settings.ircfile, 'r')
	except FileNotFoundError:
		line_prepender(settings.logfile, "CRITICAL ERROR: structure input file {0} not found!\n\n".format(settings.ircfile))
		sys.exit("CRITICAL ERROR: structure input file {0} not found!".format(settings.ircfile))		
	structure_object.seek(0)
	structure_file = (line for line in structure_object) # make generator
	settings.filetype="X"
	for line in structure_file:
		if "Gaussian, Inc." in line:
			settings.filetype="G"
			break
	if settings.filetype == "X":
		try:
			structures = structures_from_xyz(settings.ircfile)
		except StopIteration:
			line_prepender(settings.logfile, "CRITICAL ERROR: problem reading structures from file {}!\n\n".format(settings.ircfile))
			sys.exit("CRITICAL ERROR: problem reading structures from file {}!".format(settings.ircfile))	
	else:
		structures, settings = structures_from_G(settings.ircfile, settings)
#------------Get Information on Fragments/determine automatically
	#CHARGE
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator
	#check if they were extracted from file, if not set to 0 and then try to read
	try:
		settings.charge
	except AttributeError:
		settings.charge = 0
	for line in input_file:
		if len(line.split()) > 0:
			if line.split()[0].upper() == "CHARGE":
				try:
					settings.charge = int(line.split()[-1])
					break
				except IndexError:
					break
	#MULTIPLICITY
	try:
		settings.charge
	except AttributeError:
		settings.multi = 1
	for line in input_file:
		if len(line.split()) > 0:
			if line.split()[0].upper() == "MULTIPLICITY":
				try:
					settings.multi = int(line.split()[-1])
					break
				except (IndexError,ValueError) as error:
					break
	#FRAGMENT1
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	settings.frag1name = "fragment1"
	for line in input_file:
		if len(line.split()) > 1:
			if line.split()[0].upper()+line.split()[1].upper()  == "FRAGMENT1NAME":
				try:
					settings.frag1name = line.split()[-1]
					break
				except (IndexError,ValueError) as error:
					break
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	settings.frag1charge= 0
	for line in input_file:
		if len(line.split()) > 1:
			if line.split()[0].upper()+line.split()[1].upper()  == "FRAGMENT1CHARGE":
				try:
					settings.frag1charge = int(line.split()[-1])
					break
				except (IndexError,ValueError) as error:
					break						
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	settings.frag1multi = 1
	for line in input_file:
		if len(line.split()) > 1:
			if line.split()[0].upper()+line.split()[1].upper()  == "FRAGMENT1MULTIPLICITY":
				try:
					settings.frag1multi = int(line.split()[-1])
					break
				except (IndexError,ValueError) as error:
					break
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	settings.frag1energy = 0
	for line in input_file:
		if len(line.split()) > 1:
			if line.split()[0].upper()+line.split()[1].upper()  == "FRAGMENT1ENERGY":
				try:
					settings.frag1energy = float(line.split()[-1])
					break
				except (IndexError,ValueError) as error:
					break				
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	settings.frag1atoms = "auto"
	for line in input_file:
		if len(re.split("\s+|,",line)) > 2:
			if line.split()[0].upper()+line.split()[1].upper()  == "FRAGMENT1ATOMS":
				try:
					settings.frag1atoms =re.split("\s+|,",line)[3:]
					break
				except IndexError:
					break
	settings.frag1atoms = [int(x) for x in settings.frag1atoms if isInt(x)]
	if len(settings.frag1atoms) == 0:
		settings.frag1atoms  = "auto"
	#FRAGMENT2
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	settings.frag2name = "fragment2"
	for line in input_file:
		if len(line.split()) > 1:
			if line.split()[0].upper()+line.split()[1].upper()  == "FRAGMENT2NAME":
				try:
					settings.frag2name = line.split()[-1]
					break
				except IndexError:
					break
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	settings.frag2charge= 0
	for line in input_file:
		if len(line.split()) > 1:
			if line.split()[0].upper()+line.split()[1].upper()  == "FRAGMENT2CHARGE":
				try:
					settings.frag2charge = int(line.split()[-1])
					break
				except (IndexError,ValueError) as error:
					break						
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	settings.frag2multi = 1
	for line in input_file:
		if len(line.split()) > 1:
			if line.split()[0].upper()+line.split()[1].upper()  == "FRAGMENT2MULTIPLICITY":
				try:
					settings.frag2multi = int(line.split()[-1])
					break
				except (IndexError,ValueError) as error:
					break
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	settings.frag2energy = 0
	for line in input_file:
		if len(line.split()) > 1:
			if line.split()[0].upper()+line.split()[1].upper()  == "FRAGMENT2ENERGY":
				try:
					settings.frag2energy = float(line.split()[-1])
					break
				except (IndexError,ValueError) as error:
					break				
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	settings.frag2atoms = "auto"
	for line in input_file:
		if len(re.split("\s+|,",line)) > 2:
			if line.split()[0].upper()+line.split()[1].upper()  == "FRAGMENT2ATOMS":
				try:
					settings.frag2atoms =re.split("\s+|,",line)[3:]
					break
				except IndexError:
					break
	settings.frag2atoms = [int(x) for x in settings.frag2atoms if isInt(x)]
	if len(settings.frag2atoms) == 0:
		settings.frag2atoms  = "auto"
	#Determine fragment atoms if not set
	if settings.frag1atoms == "auto" and settings.frag2atoms == "auto":
		settings = auto_frag(structures,settings)
	elif settings.frag1atoms == "auto":	
		settings.frag1atoms = list(range(1,len(structures.atoms)+1))
		settings.frag1atoms = [item for item in settings.frag1atoms if item not in set(settings.frag2atoms)]
	elif settings.frag2atoms == "auto":	
		settings.frag2atoms = list(range(1,len(structures.atoms)+1))
		settings.frag2atoms = [item for item in settings.frag2atoms if item not in set(settings.frag1atoms)]	
	check_fragments(settings,structures) #checks if atom lists are coherent
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
	
#------------Get analysis information
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	settings.analysis = True
	for line in input_file:
		if len(line.split()) > 0:
			if line.split()[0].upper()  == "ANALYSIS":
				try:
					if line.split()[2].upper()  == "YES":
						settings.analysis = True
						break
					elif line.split()[2].upper()  == "NO":
						settings.analysis = False
						break
				except IndexError:
					break	
	if settings.analysis:
		settings.geo_dist = [[]]
		with open(input_filename) as input_file:
			file_contents = input_file.read() # reset generator

		for match in re.finditer(r'<distances>(.*?)</distances>', file_contents, re.IGNORECASE|re.DOTALL):
			settings.geo_dist = match.group(1).strip().split('\n')
			settings.geo_dist = [element.split() for element in settings.geo_dist]
		settings.geo_ang = []
		with open(input_filename) as input_file:
			file_contents = input_file.read() # reset generator
		for match in re.finditer(r'<angles>(.*?)</angles>', file_contents, re.IGNORECASE|re.DOTALL):
			settings.geo_ang = match.group(1).strip().split('\n')
			settings.geo_ang = [element.split() for element in settings.geo_ang]
		settings.geo_dih = []
		with open(input_filename) as input_file:
			file_contents = input_file.read() # reset generator
		for match in re.finditer(r'<dihedral>(.*?)</dihedral>', file_contents, re.IGNORECASE|re.DOTALL):
			settings.geo_dih = match.group(1).strip().split('\n')
			settings.geo_dih = [element.split() for element in settings.geo_dih]
		#Automatic determination of formed, broken bonds if requested
		if settings.geo_dist[0] == []:
			settings.geo_dist = [["auto"]]
		for element in settings.geo_dist:
			if "auto" in element:
				auto_distances = getBonds(structures)
				#replace auto in that list against the new distances
				settings.geo_dist = [ x for x in settings.geo_dist if "auto" not in x]
				settings.geo_dist.append(auto_distances[0])
				settings.geo_dist.append(auto_distances[1])
				break
		#eliminate problems with empty lists
		if settings.geo_ang[0] == []:
			settings.geo_ang = []
		if settings.geo_dih[0] == "":
			settings.geo_dih = []
		#eliminate problems with wrong inputs
		settings = check_geo(structures, settings)
#------------Get Further Setting
	#keep xyz
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	settings.keepxyz = True
	for line in input_file:
		if len(line.split()) > 2:
			if line.split()[0].upper()+line.split()[1].upper()  == "KEEPXYZ":
				try: #in case it's blank
					line.split()[3]
				except IndexError:
					settings.keepxyz = True
					break
				if line.split()[3].upper() == "YES":
					settings.keepxyz = True
				elif line.split()[3].upper() == "NO":
					settings.keepxyz = False
				else:
					settings.keepxyz = True
				break		
	#keep input files
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	settings.keepinput = True
	for line in input_file:
		if len(line.split()) > 3:
			if line.split()[0].upper()+line.split()[1].upper()+line.split()[2].upper()   == "KEEPINPUTFILES":
				try: #in case it's blank
					line.split()[4]
				except IndexError:
					settings.keepinput  = True
					break
				if line.split()[4].upper() == "YES":
					settings.keepinput = True
				elif line.split()[4].upper() == "NO":
					settings.keepinput  = False
				else:
					settings.keepinput  = True
				break				
	#keep output files
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	settings.keepoutput = True
	for line in input_file:
		if len(line.split()) > 3:
			if line.split()[0].upper()+line.split()[1].upper()+line.split()[2].upper()   == "KEEPOUTPUTFILES":
				try: #in case it's blank
					line.split()[4]
				except IndexError:
					settings.keepoutput  = True
					break
				if line.split()[4].upper() == "YES":
					settings.keepoutput = True
				elif line.split()[4].upper() == "NO":
					settings.keepoutput  = False
				else:
					settings.keepoutput  = True
				break		
	#keep log file
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	settings.keeplog = True
	for line in input_file:
		if len(line.split()) > 3:
			if line.split()[0].upper()+line.split()[1].upper()+line.split()[2].upper()   == "KEEPLOGFILE":
				try: #in case it's blank
					line.split()[4]
				except IndexError:
					settings.keeplog   = True
					break
				if line.split()[4].upper() == "YES":
					settings.keeplog  = True
				elif line.split()[4].upper() == "NO":
					settings.keeplog   = False
				else:
					settings.keeplog   = True
				break				
	#reorder
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	settings.reorder = False
	for line in input_file:
		if len(line.split()) > 2:
			if line.split()[0].upper()   == "REORDER":
				try: #in case it's blank
					line.split()[2]
				except IndexError:
					settings.reorder   = False
					break
				if line.split()[2].upper() == "YES":
					settings.reorder  = True
				elif line.split()[2].upper() == "NO":
					settings.reorder   = False
				else:
					settings.reorder   = False
				break					
	#reduce structures
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	settings.reduce = False
	for line in input_file:
		if len(line.split()) > 3:
			if line.split()[0].upper()+line.split()[1].upper()   == "REDUCESTRUCTURES":
				try: #in case it's blank
					line.split()[3]
				except IndexError:
					settings.reduce   = False
					break
				if line.split()[3].upper() == "YES":
					settings.reduce  = True
				elif line.split()[3].upper() == "NO":
					settings.reduce   = False
				else:
					settings.reduce   = False
				break	
	if settings.reduce:
		input_object.seek(0)
		input_file = (line for line in input_object) # make generator	
		settings.reduce_tresh = 0.2
		for line in input_file:
			if len(line.split()) > 3:
				if line.split()[0].upper()+line.split()[1].upper()   == "RMSDTRESHOLD":
					try: #in case it's blank
						line.split()[3]
					except IndexError:
						settings.reduce_tresh    = 0.2
						break
					try:
						settings.reduce_tresh = float(line.split()[3].upper() )
					except ValueError:
						settings.reduce_tresh = 0.2
					break
	#prepare only
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	settings.prepareonly = False
	for line in input_file:
		if len(line.split()) > 3:
			if line.split()[0].upper()+line.split()[1].upper()   == "PREPAREONLY":
				try: #in case it's blank
					line.split()[3]
				except IndexError:
					settings.prepareonly  = False
					break
				if line.split()[3].upper() == "YES":
					settings.prepareonly  = True
				elif line.split()[3].upper() == "NO":
					settings.prepareonly   = False
				else:
					settings.prepareonly   = False
				break							
	
				
#------------Get Settings for running the application
	input_object.seek(0)
	input_file = (line for line in input_object) # reset generator
	settings.input_file_extension = "com" #setting default
	for line in input_file:
		if len(line.split()) > 2: # This is needed otherwise it produces out of range error when contents of lines too short
			if line.split()[0].upper()+line.split()[1].upper()+line.split()[2].upper() == "INPUTFILEEXTENSION":
				try: #in case it's blank
					settings.input_file_extension = line.split()[4]
					break
				except IndexError:
					break
	input_object.seek(0)
	input_file = (line for line in input_object) # reset generator
	settings.output_file_extension = "log" #setting default
	for line in input_file:
		if len(line.split()) > 2: # This is needed otherwise it produces out of range error when contents of lines too short
			if line.split()[0].upper()+line.split()[1].upper()+line.split()[2].upper() == "OUTPUTFILEEXTENSION":
				try: #in case it's blank
					settings.output_file_extension = line.split()[4]
					break
				except IndexError:
					break
	#get input layout
	settings.inputlayout = ""
	with open(input_filename) as input_file:
		file_contents = input_file.read() # reset generator
	for match in re.finditer(r'<layout>(.*?)</layout>', file_contents, re.IGNORECASE|re.DOTALL):
		settings.inputlayout = match.group(1).strip()+'\n\n'
	#get run settings
	settings.submit_setting = ""
	with open(input_filename) as input_file:
		file_contents = input_file.read() # reset generator
	for match in re.finditer(r'<run_job>(.*?)</run_job>', file_contents, re.IGNORECASE|re.DOTALL):
		settings.submit_setting  = match.group(1).strip()

	return settings, structures