import sys
from os import popen

pdb_code=sys.argv[1]
pdb_code_fasta=sys.argv[2]
pdb_blast_db=sys.argv[3]
psi_blast_binary=sys.argv[4]
pdb_cif=sys.argv[5]
foldseek_pdb_lib=sys.argv[6]
foldseek_binary=sys.argv[7]
dates=sys.argv[8]
pdb_summary=sys.argv[9]
workdir=sys.argv[10]
output=sys.argv[11]

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
	
def read_fasta(fasta):
	s=""
	ID=""
	sq={}
	inf=open(fasta)
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

def prefilter(summary,code):
	res=["",""]
	inf=open(summary)
	while 1:
		line=inf.readline()
		if line=="":
			break
		line=line.split("\t")
		if line[0]==code:
			res[0]=line[1]
			res[1]=line[2].strip()
	inf.close()		
	return res	

def pdb_seq_search(sq,db,blast,tmp_):
	res={}
	log={}
	for key in sq:
		of=open(tmp_+"/"+key+".fas","w")
		of.write(">"+key+"\n"+sq[key]+"\n")	
		of.close()
		tmp=popen(blast+" -query "+tmp_+"/"+key+".fas -db "+db+" -evalue 0.0001 -max_target_seqs=50000 -num_iterations 3 -out_ascii_pssm "+tmp_+"/"+key+".pssm -out "+tmp_+"/"+key+".blast")
		content=tmp.readlines()

		try:
			log[key[:4]].append(blast+" -query "+tmp_+"/"+key+".fas -db "+db+" -evalue 0.0001 -max_target_seqs=50000 -num_iterations 3 -out_ascii_pssm "+tmp_+"/"+key+".pssm -out "+tmp_+"/"+key+".blast")		
		except Exception:
			log[key[:4]]=[]		
			log[key[:4]].append(blast+" -query "+tmp_+"/"+key+".fas -db "+db+" -evalue 0.0001 -max_target_seqs=50000 -num_iterations 3 -out_ascii_pssm "+tmp_+"/"+key+".pssm -out "+tmp_+"/"+key+".blast")	
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
			log[key[:4]].append("previous command produced an error")
	return [res,log]

def pdb_str_search(code,cif,db,foldseek,tmp_):
	res={}
	log={}
	tmp=popen(foldseek+" easy-search "+cif+" "+db+" "+tmp_+"/"+code+".foldseek "+tmp_+"/fs_tmp/ --max-seqs 50000 --exhaustive-search")
	content=tmp.readlines()
	try:
		log[code].append(foldseek+" easy-search "+cif+" "+db+" "+tmp_+"/"+code+".foldseek "+tmp_+"/fs_tmp/ --max-seqs 50000 --exhaustive-search")	
	except KeyError:
		log[code]=[]	
		log[code].append(foldseek+" easy-search "+cif+" "+db+" "+tmp_+"/"+code+".foldseek "+tmp_+"/fs_tmp/ --max-seqs 50000 --exhaustive-search")			
	try:
		inf=open(tmp_+"/"+code+".foldseek")
		while 1:
			line=inf.readline()
			if line=="":
				break
			line=line.split()

			if line[0].find("MODEL")!=-1:
				line[0]=line[0].split("_")[0]+"_"+line[0].split("_")[len(line[0].split("_"))-1]

			if float(line[2])>0.25 and int(line[7])-int(line[6])>10 and int(line[9])-int(line[8])>10:

				try:
					res[line[0]][line[1].split("-")[0]]=[line[2],int(line[7])-int(line[6]),int(line[9])-int(line[8])]
				except KeyError:
					res[line[0]]={}			
					res[line[0]][line[1].split("-")[0]]=[line[2],int(line[7])-int(line[6]),int(line[9])-int(line[8])]


	except Exception:
		log[code].append("previous command produced an error")
	return [res,log]

def selection(blast,foldseek,cut_date,summary,sq,log_b,log_f):
	res={}
	inf=open(summary)
	data={}
	while 1:
		line=inf.readline()
		if line=="":
			break
		line=line.split("\t")
		data[line[0]]=[line[1],line[2].strip()]
	inf.close()	
	
	for key in sq:
		res[key]={}
		for key2 in cut_date:
			res[key][key2]=True
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
					if method.split()[0]=="SOLUTION" or method.split()[0]=="SOLID-STATE":
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
				log_b[key[:4]].append(key+" "+key2+" "+key3+" blast "+str(blast[key][key2])+" "+str(current)+" "+date+" "+method)				
				if current==False:
					res[key][key3]=False		

							
	for key in foldseek:
		method=""
		date=""
		for key2 in foldseek[key]:
		
			try:
				method=data[key2][1]
				date=data[key2][0]
			except KeyError:
				date="1970-01-01"
				method="unknown"
			for key3 in cut_date:
				current=True
				if cut_date[key3]=="training":
					if method.split()[0]=="SOLUTION" or  method.split()[0]=="SOLID-STATE":
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
				log_f[key[:4]].append(key+" "+key2+" "+key3+" foldseek "+str(foldseek[key][key2])+" "+str(current)+" "+date+" "+method)
				if current==False:
					try:
						res[key][key3]=False
					except KeyError:
						if len(res)==1:
							for key4 in res:
								res[key4][key3]=False

	return [res,log_b,log_f]

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
		
def write_res(outp,res,log_b,log_f):
	outf=open(outp,"w")
	for key in res:
		for key2 in res[key]:
			outf.write(key+"	"+key2+"	"+str(res[key][key2])+"\n")
	outf.close()
	outf=open(outp+".log","w")
	for key in log_b:
		for i in range (0, len(log_b[key])):
			outf.write(log_b[key][i]+"\n")
	for key in log_f:
		for i in range (0, len(log_f[key])):
			outf.write(log_f[key][i]+"\n")
	outf.close()	

cutoff=read_dates(dates)
sequences=read_fasta(pdb_code_fasta)

filter_res=prefilter(pdb_summary,pdb_code)
if int(filter_res[0].split("-")[0])>=2018 and int(filter_res[0].split("-")[1])>=5:
	[blast_result,blast_log]=pdb_seq_search(sequences,pdb_blast_db,psi_blast_binary,workdir)
	[foldseek_result,foldseek_log]=pdb_str_search(pdb_code,pdb_cif,foldseek_pdb_lib,foldseek_binary,workdir)
	[search_res,blast_log,foldseek_log]=selection(blast_result,foldseek_result,cutoff,pdb_summary,sequences,blast_log,foldseek_log)
	write_res(output,search_res,blast_log,foldseek_log)
elif int(filter_res[0].split("-")[0])>2018:
	[blast_result,blast_log]=pdb_seq_search(sequences,pdb_blast_db,psi_blast_binary,workdir)
	[foldseek_result,foldseek_log]=pdb_str_search(pdb_code,pdb_cif,foldseek_pdb_lib,foldseek_binary,workdir)
	[search_res,blast_log,foldseek_log]=selection(blast_result,foldseek_result,cutoff,pdb_summary,sequences,blast_log,foldseek_log)
	write_res(output,search_res,blast_log,foldseek_log)	
elif filter_res[1].split()[0]=="SOLUTION" or filter_res[1].split()[0]=="SOLID-STATE":
	[blast_result,blast_log]=pdb_seq_search(sequences,pdb_blast_db,psi_blast_binary,workdir)
	[foldseek_result,foldseek_log]=pdb_str_search(pdb_code,pdb_cif,foldseek_pdb_lib,foldseek_binary,workdir)
	[search_res,blast_log,foldseek_log]=selection(blast_result,foldseek_result,cutoff,pdb_summary,sequences,blast_log,foldseek_log)
	search_res=correct(search_res,cutoff,filter_res[0])
	write_res(output,search_res,blast_log,foldseek_log)
else:
	search_res=filler(cutoff,sequences)
	write_res(output,search_res,{'-':[]},{'-':[]})


