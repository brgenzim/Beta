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
			ID=l.split("|")[1]
			s=""
		else:
			s=s+l.strip()
	
	inf.close()
	return sq

parameters=read_param(parameters_file)
seqs=readfasta(parameters["SPROT_FASTA"])


outf=open(parameters["SEQ_SEARCH_RESULTS"],"w")
for key in seqs:

	
	inf=open(parameters["WORKDIR"]+"/"+'up_'+key[1:3]+'/'+key+'.result')
	while 1:
		line=inf.readline()
		if line=="":
			break
		line=line.split()
		
		if line[2].strip()=="True":

			outf.write(line[0]+"\t"+line[1]+"\n")


													
														
	inf.close()
	
outf.close()	
