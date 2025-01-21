
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


uniprot_results={}
interactions={}

inf=open(parameters["SEQ_SEARCH_RESULTS"])
while 1:
	line=inf.readline()
	if line=="":
		break
	line=line.split("\t")
	try:
		uniprot_results[line[0]][line[1]]=""
	except KeyError:
		uniprot_results[line[0]]={}
		uniprot_results[line[0]][line[1]]=""			
inf.close()


inf=open(parameters["BIOGRID"])
while 1:
	line=inf.readline()
	if line=="":
		break
	line=line.split("\t")
	
	ids=line[2].split("|")
	AC1=""
	AC2=""
	for i in range (0, len(ids)):
		if ids[i].find("swiss-prot")!=-1:
			AC1=ids[i].split(":")[1]

	ids=line[3].split("|")
	for i in range (0, len(ids)):
		if ids[i].find("swiss-prot")!=-1:
			AC2=ids[i].split(":")[1]				
	
	try:
		pr=uniprot_results[AC1]
		pr=uniprot_results[AC2]		
		try:
			pr=interactions[AC1]
			try:
			
				interactions[AC1][AC2]=""
			except KeyError:
				interactions[AC1]={}
				interactions[AC1][AC2]=""		
		except KeyError:
			try:
			
				interactions[AC2][AC1]=""
			except KeyError:
				interactions[AC2]={}
				interactions[AC2][AC1]=""		
	except KeyError:
		pass
inf.close()



outf=open(parameters["SEQ_SEARCH_INTERACTIONS_RESULT"],"w")
c=0
for key in interactions:
	for key2 in interactions[key]:
		for key3 in uniprot_results[key]:
			for key4 in uniprot_results[key2]:		
				if key3==key4:
					outf.write(key+"	"+key2+"	"+key3)
outf.close()

