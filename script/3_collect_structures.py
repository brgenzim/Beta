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
		if line.find("=")!=-1:
			line=line.split("=")
			p[line[0]]=line[1].strip()
	return p
	
	
def readsifts(x):
	pdb_={}
	inf=open(x)
	line=inf.readline()
	line=inf.readline()
	while 1:
		line=inf.readline()
		if line=="":
			break
		try:
			pdb_[line.split()[0]].append(line.split()[1])
		except KeyError:
			pdb_[line.split()[0]]=[]
			pdb_[line.split()[0]].append(line.split()[1])			
	inf.close()
	return pdb_

def readfasta(x):
	s=""
	ID=""
	sq={}
	inf=open(x)
	while 1:
		l=inf.readline()
		if l=="":
			sq[ID]=s
			break
		if l[0]==">":
			if ID!="":
				sq[ID]=s
			ID=l[1:].split()[0]
			s=""
		else:
			s=s+l.strip()
	
	inf.close()
	return sq
def prefilter(summary):
	res={}
	res2={}	
	res3={}		
	inf=open(summary)
	while 1:
		line=inf.readline()
		if line=="":
			break
		line=line.split("\t")
		res[line[0]]=line[1]
		res2[line[0]]=line[2]
		res3[line[0]]=line[3].strip().replace("?","0")		
		
	inf.close()		
	return [res,res2,res3]
def read_dates(cut_date):
	dt={}
	inf=open(cut_date)
	while 1:
		line=inf.readline()
		if line=="":
			break
		line=line.split()
		dt[line[0]]=line[1].strip()
	return dt
		
parameters=read_param(parameters_file)	
pdb_ids=readsifts(parameters["SIFTS"])
seqs=readfasta(parameters["PDB_FAS"])
[release,method,resolution]=prefilter(parameters["WORKDIR"]+"/pdb_summary.txt")
dates=read_dates(parameters["CUTOFF"])
max_resolution=float(parameters["RESOLUTION"])
outf=open(parameters["STR_SEARCH_RESULTS"],"w")
for key in pdb_ids:
	inf=open(parameters["WORKDIR"]+'pdb_'+key[1:3]+'/'+key+'.result')
	while 1:
		line=inf.readline()
		if line=="":
			break
		line=line.split()
		
		if line[2].strip()=="True":

			
			if method[line[0][:4]].find("NMR")==-1:			
				
				if int(release[line[0][:4]].split("-")[0])>int(line[1].split("-")[0]) and float(resolution[line[0][:4]])<=max_resolution:
					outf.write(line[0]+"\t"+line[1]+"\n")				

				elif int(release[line[0][:4]].split("-")[0])==int(line[1].split("-")[0]) and int(release[line[0][:4]].split("-")[1])>int(line[1].split("-")[1]) and float(resolution[line[0][:4]])<=max_resolution:
					outf.write(line[0]+"\t"+line[1]+"\n")


			else:
				if dates[line[1]]=="template":
					if int(release[line[0][:4]].split("-")[0])>int(line[1].split("-")[0]):
						outf.write(line[0]+"\t"+line[1]+"\n")				
								
					elif int(release[line[0][:4]].split("-")[0])==int(line[1].split("-")[0]) and int(release[line[0][:4]].split("-")[1])>int(line[1].split("-")[1]):
						outf.write(line[0]+"\t"+line[1]+"\n")					
				else:
					outf.write(line[0]+"\t"+line[1]+"\n")										
														
	inf.close()

outf.close()	
