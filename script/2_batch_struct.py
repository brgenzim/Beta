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
	
	
parameters=read_param(parameters_file)
pdb_ids=readsifts(parameters["SIFTS"])
seqs=readfasta(parameters["PDB_FAS"])

start=time.time()
pdb_count=0
used={}
used2={}

for key in pdb_ids:
	pdb_count+=1	
	progresser="0"
	if pdb_count%100==0:
		while 1:
			tmp=popen("squeue -u "+parameters["USER"])
			content=tmp.readlines()
			if len(content)<200:
				break
			time.sleep(10)
			if str((pdb_count/len(pdb_ids))*100)[:5]!=progresser:
				print("progress: "+str((pdb_count/len(pdb_ids))*100)[:5]+'%')
				progresser=str((pdb_count/len(pdb_ids))*100)[:5]
			
	try:
		pr=used[key[1:3]]
	except KeyError:
		tmp=popen("mkdir "+parameters["WORKDIR"]+"pdb_"+key[1:3])
		content=tmp.readlines()
		used[key[1:3]]=""

		
	outf=open(parameters["WORKDIR"]+"pdb_"+key[1:3]+"/"+key+".fas","w")
	for i in range (0, len(pdb_ids[key])):
		try:
			outf.write(">"+key+"_"+pdb_ids[key][i]+"\n")
			outf.write(seqs[key+"_"+pdb_ids[key][i]]+"\n")		
		except KeyError:
			pass
	outf.close()
	try:
		pr=used2[key]
	except KeyError:
		tmp=popen("mkdir "+parameters["WORKDIR"]+"pdb_"+key[1:3]+"/"+key)
		content=tmp.readlines()	
		used2[key]=""
	tmp=popen('zless '+parameters["PDB"]+"/pdb/structures/divided/mmCIF/"+key[1:3]+'/'+key+'.cif.gz | grep "" >'+parameters["WORKDIR"]+'pdb_'+key[1:3]+'/'+key+'/'+key+'.cif')	
	content=tmp.readlines()	
	


	tmp=popen('sbatch --wrap="python3 '+parameters["SCRIPT1"]+' '+key+' '+parameters["WORKDIR"]+'pdb_'+key[1:3]+'/'+key+'.fas '+parameters["PDB_FAS"]+' '+parameters["PSIBLAST"]+' '+parameters["WORKDIR"]+'pdb_'+key[1:3]+'/'+key+'/'+key+'.cif'+' '+parameters['PDB_FOLDSEEK']+' foldseek '+parameters["CUTOFF"]+' '+parameters["WORKDIR"]+'/pdb_summary.txt '+parameters["WORKDIR"]+' '+parameters["WORKDIR"]+'pdb_'+key[1:3]+'/'+key+'.result" --mem 4G -c 2 -o '+parameters["WORKDIR"]+'pdb_'+key[1:3]+'/'+key+'.slurmlog -e '+parameters["WORKDIR"]+'pdb_'+key[1:3]+'/'+key+'.slurmerror')
	content=tmp.readlines()
	
	

	

