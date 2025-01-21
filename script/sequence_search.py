import sys
from os import popen


uniprot_ac=sys.argv[1]
uniprot_fasta=sys.argv[2]
pdb_blast_db=sys.argv[3]
psi_blast_binary=sys.argv[4]
sifts=sys.argv[5]
dates=sys.argv[6]
pdb_result=sys.argv[7]
pdb_summary=sys.argv[8]
workdir=sys.argv[9]
output=sys.argv[10]

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
	
def read_fasta(fasta,ac):
	s=""
	ID=""
	sq={}
	inf=open(fasta)
	while 1:
		l=inf.readline()
		if l=="":
			sq[ac]=s
			break
		if l[0]==">":
			if ID!="":
				sq[ac]=s
			ID=l[1:].strip()
			s=""
		else:
			s=s+l.strip()
	
	inf.close()
	return sq

def readsifts(ac,x):
	pdb_=[]
	inf=open(x)
	line=inf.readline()
	line=inf.readline()
	while 1:
		line=inf.readline()
		if line=="":
			break
		line=line.split()
		if line[2]==ac:
			pdb_.append(line[0]+"_"+line[1])

	inf.close()
	pdb_=list(set(pdb_))
	return pdb_

def read_pdb_results(prec):
	res={}
	inf=open(prec)
	while 1:
		line=inf.readline()
		if line=="":
			break
		line=line.split()
		try:
			res[line[0]][line[1]]=""
		except KeyError:
			res[line[0]]={}		
			res[line[0]][line[1]]=""		
	inf.close()	
	return res
def prefilter(res,precalc,code,ac):

	for key in res:
		for key2 in res[key]:
			try:

				pr=precalc[code][key2]

			except KeyError:
				print(code+" "+key2)	

				res[key][key2]=False
	return res
def prefilter2(res):
	check=False	
	for key in res:
		for key2 in res[key]:
			if res[key][key2]==True:
				check=True
	return check
		
def seq_search(sq,db,blast,tmp_):
	res={}
	log={}
	for key in sq:
		of=open(tmp_+"/"+key+".fas","w")
		of.write(">"+key+"\n"+sq[key]+"\n")	
		of.close()
		tmp=popen(blast+" -query "+tmp_+"/"+key+".fas -db "+db+" -evalue 0.0001 -max_target_seqs=50000 -num_iterations 3 -out_ascii_pssm "+tmp_+"/"+key+".pssm -out "+tmp_+"/"+key+".blast")
		content=tmp.readlines()

		try:
			log[key].append(blast+" -query "+tmp_+"/"+key+".fas -db "+db+" -evalue 0.0001 -max_target_seqs=50000 -num_iterations 3 -out_ascii_pssm "+tmp_+"/"+key+".pssm -out "+tmp_+"/"+key+".blast")		
		except Exception:
			log[key]=[]		
			log[key].append(blast+" -query "+tmp_+"/"+key+".fas -db "+db+" -evalue 0.0001 -max_target_seqs=50000 -num_iterations 3 -out_ascii_pssm "+tmp_+"/"+key+".pssm -out "+tmp_+"/"+key+".blast")	
		inf=open(tmp_+"/"+key+".blast")						

		try:
			while 1:
				line=inf.readline()
				if line=="":
					break
				if line[0]==">":

					ID=line[1:5]
					line=inf.readline()
					line=inf.readline()
					line=inf.readline()
					line=inf.readline()
					if len(line.split())==0:
						while 1:
							line=inf.readline()
							try:
								if line.split()[0]=="Identities":
									break					
							except IndexError:
								pass											
					if line.split()[0]!="Identities":
						while 1:
							line=inf.readline()
							try:
								if line.split()[0]=="Identities":
									break					
							except IndexError:
								pass											
					line=line.split()



					if int(line[3].replace("(","").replace("%)","").replace(",",""))>20 and int(line[2].split("/")[1])>10:
						try:
							res[key][ID]=[int(line[3].replace("(","").replace("%)","").replace(",","")),int(line[2].split("/")[1])]		
						except KeyError:
							res[key]={}
							res[key][ID]=[int(line[3].replace("(","").replace("%)","").replace(",","")),int(line[2].split("/")[1])]							

		except Exception:
			log[key].append("previous command produced an error")
	return [res,log]



def selection(blast,cut_date,summary,sq,log_b,res):

	inf=open(summary)
	data={}
	while 1:
		line=inf.readline()
		if line=="":
			break
		line=line.split("\t")
		data[line[0]]=[line[1],line[2].strip()]
	inf.close()	
	

	for key in blast:

		method=""
		date=""
		for key2 in blast[key]:
			try:
				method=data[key2][1]
				date=data[key2][0]
			except KeyError:
				date="1970-01-01"
				method="unknown"
			for key3 in cut_date:
				current=True
				if cut_date[key3]=="training":
					if method.split()[0]=="SOLUTION":
						pass
					else:
						if int(date.split("-")[0])<=int(key3.split("-")[0]) and int(date.split("-")[1])<=int(key3.split("-")[1]):
							current=False
						elif int(date.split("-")[0])<int(key3.split("-")[0]):					
							current=False
				else:

					if int(date.split("-")[0])<=int(key3.split("-")[0]) and int(date.split("-")[1])<=int(key3.split("-")[1]):
						current=False
					elif int(date.split("-")[0])<int(key3.split("-")[0]):					
						current=False	
				log_b[key].append(key+" "+key2+" "+key3+" blast "+str(blast[key][key2])+" "+str(current)+" "+date+" "+method)				
				if current==False:
					res[key][key3]=False		


	
	return [res,log_b]

def filler(cut_date,sq):
	res={}
	for key in sq:
		res[key]={}
		for key2 in cut_date:
			res[key][key2]=False
	return res
	
def correct(res,cut_date,pdb_date):

	for key in res:
		for key2 in res[key]:
			if cut_date[key2]=="template":
		
				if int(pdb_date.split("-")[0])<=int(key2.split("-")[0]) and int(pdb_date.split("-")[1])<=int(key2.split("-")[1]):
					res[key][key2]=False
				elif int(pdb_date.split("-")[0])<int(key2.split("-")[0]):					
					res[key][key2]=False					

	
	return res
		
def write_res(outp,res,log_b):
	outf=open(outp,"w")
	for key in res:
		for key2 in res[key]:
			outf.write(key+"	"+key2+"	"+str(res[key][key2])+"\n")
	outf.close()
	outf=open(outp+".log","w")
	for key in log_b:
		for i in range (0, len(log_b[key])):
			outf.write(log_b[key][i]+"\n")

	outf.close()	

pdb_code=readsifts(uniprot_ac,sifts)
cutoff=read_dates(dates)
sequences=read_fasta(uniprot_fasta,uniprot_ac)
precalculated=read_pdb_results(pdb_result)
search_res={}
search_res[uniprot_ac]={}
for key in cutoff:
	search_res[uniprot_ac][key]=True
print(search_res)
for i in range (0, len(pdb_code)):
	search_res=prefilter(search_res,precalculated,pdb_code[i],uniprot_ac)
print(search_res)	
if prefilter2(search_res)==True:
	[blast_result,blast_log]=seq_search(sequences,pdb_blast_db,psi_blast_binary,workdir)
	print(search_res)		
	[search_res,blast_log]=selection(blast_result,cutoff,pdb_summary,sequences,blast_log,search_res)
	print(search_res)		
	write_res(output,search_res,blast_log)
	
else:
	write_res(output,search_res,{'-':[]})
