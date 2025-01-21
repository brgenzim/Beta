parameters_file="../parameters.txt"

def read_param(x):
	p={}
	inf=open(x)
	while 1:
		line=inf.readline()
		if line=="":
			break
		if line.find("=")!=-1:
			line=line.split("=")
			p[line[0]]=line[1].strip()
	return p
	
parameters=read_param(parameters_file)
	
pdbs={}
chains={}
inf=open(parameters["STR_SEARCH_RESULTS"])
while 1:
	line=inf.readline()
	if line=="":
		break
	line=line.split()
	pdbs[line[0][:4]]=""
	try:
		chains[line[0][:4]][line[0].split("_")[1]][line[1].strip()]=""
	except KeyError:
		try:
			chains[line[0][:4]][line[0].split("_")[1]]={}
			chains[line[0][:4]][line[0].split("_")[1]][line[1].strip()]=""
		except KeyError:
			chains[line[0][:4]]={}
			chains[line[0][:4]][line[0].split("_")[1]]={}
			chains[line[0][:4]][line[0].split("_")[1]][line[1].strip()]=""				
inf.close()

	

outf=open(parameters["STR_SEARCH_INTERACTIONS_RESULT"],"w")
for key in pdbs:
	
	inf=open(parameters["WORKDIR"]+"pdb_"+key[1:3]+"/"+key+".voronota")
	pairs=[]
	while 1:
		line=inf.readline()
		if line=="":
			break
		pairs.append([line.split()[0],line.split()[1].strip()])
	inf.close()



	for key2 in chains[key]:
		for i in range (0, len(pairs)):
			dates1={}
			dates2={}			

			try:

				for key3 in chains[key][pairs[i][0]]:
					dates1[key3]=""

			except KeyError:
				pass
			try:
				for key3 in chains[key][pairs[i][1]]:
					dates2[key3]=""

			except KeyError:
				pass				

			for key3 in dates1:
				try:
					pr=dates2[key3]
					outf.write(key+" "+pairs[i][0]+" "+pairs[i][1]+" "+key3+"\n")
					

				except KeyError:
					pass

outf.close()
