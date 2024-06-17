from os import popen
import sys

PDB=sys.argv[1]
INF=sys.argv[2]
OUTF=sys.argv[3]


tmp=popen('zless '+INF+'/'+PDB+'-assembly1.cif.gz | grep ""')
content=tmp.readlines()


interactions={}
if len(content)<50000:
	tmp=popen('voronota-contacts -i '+INF+'/'+PDB+'-assembly1.cif.gz')
	content4=tmp.readlines()
	for i in range (0, len(content4)):
		if content4[i].find("solvent")==-1:
			chain1=""
			chain2=""
			for j in range (0, len(content4[i])-2):
				if content4[i][j]=="c" and content4[i][j+1]=="<":
					if chain1=="":	
						posi=j+2
						while 1:
							if content4[i][posi]==">":							
								chain1=content4[i][j+2:posi]	
								break
							posi+=1					
					else:
						posi=j+2
						while 1:
							if content4[i][posi]==">":							
								chain2=content4[i][j+2:posi]	
								break
							posi+=1	

			if chain1!=chain2:
		
				try:
					pr=interactions[chain2][chain1]
				except KeyError:
					try:
						pr=interactions[chain1][chain2]
					except KeyError:					
						try:
							interactions[chain1][chain2]=""
						except KeyError:
							interactions[chain1]={}						
							interactions[chain1][chain2]=""
						
							
outf=open(OUTF+"/"+PDB+".voronota","w")
for key in interactions:
	for key2 in interactions[key]:
		outf.write(key+"	"+key2+"\n")
outf.close()		
