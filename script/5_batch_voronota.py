from os import popen
import time

parameters_file="../parameters.txt"

def read_param(x):
	p={}
	inf=open(x)
	while 1:
		line=inf.readline()
		if line=="":
			break
		line=line.split("#")[0]
		line=line.split("=")
		p[line[0]]=line[1]
	return p
	
parameters=read_param(parameters_file)

pdb={}
inf=open(parameters["STR_SEARCH_RESULTS"])
while 1:
	line=inf.readline()
	if line=="":
		break
	pdb[line[:4]]=""
inf.close()
pdb_count=1
for key in pdb:
	if pdb_count%100==0:
		while 1:
			tmp=popen("squeue -u"+USER)
			content=tmp.readlines()
			if len(content)<200:
				break
			time.sleep(10)		

		tmp=popen('sbatch --wrap="python3 voronota.py '+key+' '+parameters["ASSEMBLIES"]+' '+parameters["WORKDIR"]+"/pdb_"+key[1:3]+'/" -o '+parameters["WORKDIR"]+"/pdb_"+key[1:3]+'/'+key+'.voronota.slurmlog -e '+parameters["WORKDIR"]+"/pdb_"+key[1:3]+'/'+key+'.voronota.slurmerror')
		content=tmp.readlines()

	pdb_count+=1

