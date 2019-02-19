#fragments test
class settings:
	pass
class structures:
	pass

import numpy as np


periodic_table = ["","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
    "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl",
    "Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"]
	
#Covalent radii taken from DOI: 10.1039/b801115j
#Everything beyond Cm was set to 1.80
covalent_radii = [0.00,0.32,0.28,1.28,0.96,0.84,0.76,0.71,0.66,0.57,0.58,1.66,1.41,1.21,1.11,1.07,1.05,1.02,1.06,2.03,1.76,1.70,1.60,1.53,1.39,1.61,1.52,1.50,1.24,1.32,1.22,1.22,1.20,1.19,1.20,1.20,1.16,2.20,1.95,1.90,1.75,
    1.64,1.54,1.47,1.46,1.42,1.39,1.45,1.44,1.42,1.39,1.39,1.38,1.39,1.40,2.44,2.15,2.07,2.04,2.03,2.01,1.99,1.98,1.98,1.96,1.94,1.92,1.92,1.89,1.90,1.87,1.87,1.75,1.70,1.62,1.51,1.44,1.41,1.36,1.36,1.32,1.45,
    1.46,1.48,1.40,1.50,1.50,2.60,2.21,2.15,206,2.00,1.96,1.90,1.87,180,169,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80]
	
def isInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False	
#takes xyz coordinates as numpy array and returns distance
def getDistance(atom1, atom2):
	distance = 0.00
	distance = np.sqrt((atom1[0]-atom2[0])**2+(atom1[1]-atom2[1])**2+(atom1[2]-atom2[2])**2)
	return distance

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
	#then check again and so on until nothing else changes, iterative process. Take length of fragments for that since it converged if amount of fragments doesn't changes
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
	
def auto_frag(structures,settings):
	list_of_fragments = []
	for i in range(0, len(structures.xyz)):
		list_of_fragments =  list_of_fragments + [getFragments(structures.xyz[i], structures.atoms)]
	n_fragments = []
	print(list_of_fragments)
	for i in range(0,len(list_of_fragments)):
		n_fragments = n_fragments + [len(list_of_fragments[i])]
	#the next steps only work in cases where two fragments were found
	try:
		indices = [i for i, x in enumerate(n_fragments) if x == 2]
		counts = [list_of_fragments.count(list_of_fragments[x]) for x in indices]
		fragments = list_of_fragments[indices[counts.index(max(counts))]]
		settings.frag1atoms = np.sort(fragments[0])
		settings.frag2atoms = np.sort(fragments[1])
	except:
		pass
	return settings

settings.ircfile = "test.xyz"
structures = structures_from_xyz(settings.ircfile)
settings = auto_frag(structures, settings)
print(settings.frag1atoms)
print(settings.frag2atoms)