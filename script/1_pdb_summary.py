from os import popen
import gemmi
from gemmi import cif

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

def readsifts(sift):
	pdb_={}
	inf=open(sift)
	line=inf.readline()
	line=inf.readline()
	while 1:
		line=inf.readline()
		if line=="":
			break
		pdb_[line.split()[0]]=""
	inf.close()
	return pdb_

def out(outp,pdb_,prefix):
	outf=open(outp,"w")
	for key in pdb_:
		doc = cif.read(prefix+'/pdb/structures/divided/mmCIF/'+key[1:3]+'/'+key+'.cif.gz')
		block = doc.sole_block()
		method="unknown"
		date="unknown"
		try:
			method=block.find_value("_exptl.method")
			method=method.replace("'","")			
		except Exception:
			pass
		try:
			date=block.find_value("_pdbx_database_status.recvd_initial_deposition_date")
		except Exception:
			pass	
		if str(date)=="None":
			date="unknown"
		if str(method)=="None":
			method="unknown"			
		outf.write(key+"\t"+date+"\t"+method+"\n")
	outf.close()



parameters=read_param(parameters_file)
pdb_ids=readsifts(parameters["SIFTS"])
out(parameters["PDB_SUMMARY_FILE"],parameters["PDB"])
