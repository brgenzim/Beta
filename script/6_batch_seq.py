from os import popen
import time



parameters_file="../parameters.txt"

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
			ID=l[1:].split("|")[1]
			s=""
		else:
			s=s+l.strip()
	
	inf.close()
	return sq

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

seqs=readfasta(parameters["SPROT_FASTA"])

start=time.time()
count=0
used={}
progresser="0"
for key in seqs:
	count+=1
	if count%100==0:
		while 1:
			tmp=popen("squeue -u "+parameters["USER"])
			content=tmp.readlines()
			if len(content)<200:
				break
			time.sleep(10)
			if str((count/len(seqs))*100)[:5]!=progresser:
				print("progress: "+str((count/len(seqs))*100)[:5]+'%')
				progresser=str((count/len(seqs))*100)[:5]				
	try:
		pr=used[key[1:3]]
	except KeyError:
		tmp=popen("mkdir "+parameters["WORKDIR"]+"up_"+key[1:3])
		content=tmp.readlines()
		used[key[1:3]]=""
		
	outf=open(parameters["WORKDIR"]+"up_"+key[1:3]+"/"+key+".fas","w")
	outf.write(">"+key+"\n")
	outf.write(seqs[key]+"\n")		
	outf.close()
	tmp=popen('sbatch --wrap="python3 '+parameters["SCRIPT2"]+' '+key+' '+parameters["WORKDIR"]+'up_'+key[1:3]+'/'+key+'.fas '+parameters["PDB_FAS"]+' '+parameters["PSIBLAST"]+' '+parameters["SIFTS"]+' '+parameters["CUTOFF"]+' '+parameters["STR_SEARCH_RESULTS"]+' '+parameters["WORKDIR"]+'/pdb_summary.txt '+parameters["WORKDIR"]+'up_'+key[1:3]+'/ '+parameters["WORKDIR"]+'up_'+key[1:3]+'/'+key+'.result" --mem 4G -o '+parameters["WORKDIR"]+'up_'+key[1:3]+'/'+key+'.slurmlog -e '+parameters["WORKDIR"]+'up_'+key[1:3]+'/'+key+'.slurmerror')
	content=tmp.readlines()


		
