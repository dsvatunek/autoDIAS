class settings:
	pass
class structures:
	pass
	
periodic_table = ["","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
    "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl",
    "Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"]

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
	strcutures.xyz.append(XYZ) #append first structure to structures list
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


def parse_in(input_filename):
	import re
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
		sys.exit()	
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
	settings.logfile= settings.name + "_log.txt"
#------------Parse Structures
	#find out what kind of file we are dealing with. See if it's a gaussian file
	structure_object = open(settings.ircfile, 'r')
	structure_object.seek(0)
	structure_file = (line for line in input_object) # make generator
	settings.filetype="X"
	for line in structure_file:
		if "Gaussian, Inc." in line:
			settings.filetype="G"
	if setting.filetype == "X"
		structures = structures_from_xyz(settings.ircfile)
	else:
		structures = structures_from_G(settings.ircfile)	
#------------Get Information on Fragments/determine automatically
	#CHARGE
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator
	#check if they were extracted from file, if not set to 0 and then try to read
	try:
		settings.charge
	except NameError:
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
	except NameError:
		settings.multi = 1
	for line in input_file:
		if len(line.split()) > 0:
			if line.split()[0].upper() == "MULTIPLICITY":
				try:
					settings.multi = int(line.split()[-1])
					break
				except IndexError:
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
				except IndexError:
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
				except IndexError:
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
				except IndexError:
					break
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	settings.frag1energy = 0
	for line in input_file:
		if len(line.split()) > 1:
			if line.split()[0].upper()+line.split()[1].upper()  == "FRAGMENT1ENERGY":
				try:
					settings.frag1multi = float(line.split()[-1])
					break
				except IndexError:
					break				
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	settings.frag1atoms = "auto"
	for line in input_file:
		if len(re.split("\s+|,",line)) > 1:
			if line.split()[0].upper()+line.split()[1].upper()  == "FRAGMENT1ATOMS":
				try:
					settings.frag1atoms =re.split("\s+|,",line)[3:]
					break
				except IndexError:
					break				
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
				except IndexError:
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
				except IndexError:
					break
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	settings.frag2energy = 0
	for line in input_file:
		if len(line.split()) > 1:
			if line.split()[0].upper()+line.split()[1].upper()  == "FRAGMENT2ENERGY":
				try:
					settings.frag2multi = float(line.split()[-1])
					break
				except IndexError:
					break				
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	settings.frag2atoms = "auto"
	for line in input_file:
		if len(re.split("\s+|,",line)) > 1:
			if line.split()[0].upper()+line.split()[1].upper()  == "FRAGMENT2ATOMS":
				try:
					settings.frag2atoms =re.split("\s+|,",line)[3:]
					break
				except IndexError:
					break
	#Determine fragment atoms if not set
	
#------------Get analysis information
	input_object.seek(0)
	input_file = (line for line in input_object) # make generator	
	settings.analysis = True
	for line in input_file:
		if len(line.split()) > 0:
			if line.split()[0].upper()  == "ANALYSIS":
				try:
					if settings.analysis.upper() == "YES":
						settings.analysis = True
					elif settings.analysis.upper() == "NO":
						settings.analysis = False
				except IndexError:
					break	
	if settings.analysis:
		settings.geo_dist = []
		with open(input_filename) as input_file:
			file_contents = input_file.read() # reset generator
		try:
			for match in re.finditer(r'<distances>(.*?)</distances>', file_contents, re.IGNORECASE|re.DOTALL):
				settings.geo_dist = match.group(1).strip().split('\n')
				settings.geo_dist = [element.split() for element in settings.geo_dist]
		settings.geo_ang = []
		with open(input_filename) as input_file:
			file_contents = input_file.read() # reset generator
		try:
			for match in re.finditer(r'<angles>(.*?)</angles>', file_contents, re.IGNORECASE|re.DOTALL):
				settings.geo_ang = match.group(1).strip().split('\n')
				settings.geo_ang = [element.split() for element in settings.geo_ang]
		settings.geo_ang = []
		with open(input_filename) as input_file:
			file_contents = input_file.read() # reset generator
		try:
			for match in re.finditer(r'<dihedral>(.*?)</dihedral>', file_contents, re.IGNORECASE|re.DOTALL):
				settings.geo_dih = match.group(1).strip().split('\n')
				settings.geo_dih = [element.split() for element in settings.geo_dih]
		#Automatic determination of formed, broken bonds if requested
		for element in settings.geo_dist:
			if "auto" in element:
				auto_distances = auto_dist(structures)
				#replace auto in that list against the new distances
				settings.geo_dist = [ x for x in settings.geo_dist if "auto" not in x]
				settings.geo_dist.append(auto_distances[0])
				settings.geo_dist.append(auto_distances[1])
				break
#------------Get Further Setting

#------------Get Settings for running the application




	
	return settings, structures