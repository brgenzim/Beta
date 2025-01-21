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
		if line.find("=")!=-1:
			line=line.split("=")
			p[line[0]]=line[1].strip()
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
	tmp=popen("mkdir "+parameters["WORKDIR"])
	content=tmp.readlines()
	outf=open(outp+"/pdb_summary.txt","w")
	pdb_count=0
	progresser="0"
	for key in pdb_:
		pdb_count+=1
		if pdb_count%1000==0:
			if str((pdb_count/len(pdb_))*100)[:5]!=progresser:
				print("progress: "+str((pdb_count/len(pdb_))*100)[:5]+'%')
				progresser=str((pdb_count/len(pdb_))*100)[:5]
		doc = cif.read(prefix+'/pdb/structures/divided/mmCIF/'+key[1:3]+'/'+key+'.cif.gz')
		block = doc.sole_block()
		method="unknown"
		date="unknown"
		resolution="unknown"
		method=block.find_value("_exptl.method")	
		date=block.find_value("_pdbx_database_status.recvd_initial_deposition_date")
		resolution=block.find_value("_refine.ls_d_res_high")
		if str(date)=="None" or str(date)=="NoneType":
			date="unknown"
		if str(method)=="None" or str(method)=="NoneType":
			try:
				method=block.find_loop("_exptl.method")[0]
				if str(method)=="None" or str(method)=="NoneType":
					method="unknown"
			except IndexError:
				method="unknown"			
		if str(resolution)=="None" or str(resolution)=="NoneType":	
			try:
				resolution=block.find_loop("_refine.ls_d_res_high")[0]
				if str(resolution)=="None" or str(resolution)=="NoneType":
					resolution="unknown"			
			except IndexError:
				resolution="unknown"	
		method=method.replace("'","")						
		resolution=resolution.replace("'","")
		if str(resolution)=="unknown" or str(method)=="ELECTRON MICROSCOPY":	
			resolution=block.find_value("_em_3d_reconstruction.resolution")
			if str(resolution)=="None" or str(resolution)=="NoneType":
				try:
					resolution=block.find_loop("_em_3d_reconstruction.resolution")[0]
					if str(resolution)=="None" or str(resolution)=="NoneType":
						resolution="unknown"			
				except IndexError:
					resolution="unknown"								
		resolution=resolution.replace("'","")					
		outf.write(key+"\t"+date+"\t"+method+"\t"+resolution+"\n")
	outf.close()

parameters=read_param(parameters_file)
pdb_ids=readsifts(parameters["SIFTS"])
out(parameters["WORKDIR"],pdb_ids,parameters["PDB"])
