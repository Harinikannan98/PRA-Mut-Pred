#!/opt/websites/anaconda/envs/hari_37/bin/python
'''#!/usr/bin/env python3.7'''
import cgi, os
from functools import total_ordering
import sys
import time
import shutil
import subprocess
import math
import numpy as np
from os import path
import ssl
#ssl._create_default_https_context = ssl._create_unverified_context
import uuid
import cgitb;
cgitb.enable()
import timeit
start = timeit.default_timer()
import shutil
#import wget
import glob
import Bio.PDB as bpdb
from Bio.PDB import is_aa
from Bio.PDB import PDBParser, PDBIO, Select
import urllib
import os
import numpy as np
import re
import pandas as pd
import math
from Bio import PDB
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning
'''
form = cgi.FieldStorage()
pdb_id=form.getvalue('fname')
mutant=form.getvalue('input_mutations')
if form['pdbf'].filename:          
	fileitem=form['pdbf'].filename
else:
	fileitem=b""
if form['mutf'].filename:          
	mutfile=form['mutf'].filename
else:
	mutfile=b""
pdb_id=pdb_id.lower()
'''
do=input("Enter the type of input: 'pdb-id' or 'pdb-file'")
if do=="pdb-id":
	inn=input("Enter the PDB-ID")
	pdb_id=inn
	#chain=""
	#dchain=""
	pdb_file=b""
	pdb_id=pdb_id.lower()
if do=='pdb-file':
	inn=input("Enter the PDB file name in current folder")
	pdb_file=inn
	pdb_id=""
	#chain=""
	#dchain=""
mutant=input("Please the mutations in the following format 'chain:wt-position-mutation'.(For eg. A:F51A For mutations has to be made in multiple chains of homomer, use as A-B:H49D)")
mutfile=b""
'''
form = cgi.FieldStorage()
pdb_id=form.getvalue('fname')
chain=form.getvalue('chain')
dchain=form.getvalue('lchain')

if form['pdbf'].filename:          
	fileitem=form['pdbf'].filename
else:
'''
#fileitem=b""
#rna_strand=form.getvalue('rna_strand')
pdb_id=pdb_id.lower()
#import sys
sys.stdout.flush()
print ('Content-Type: text/html\r\n')
#print('/r/n')
#os.system("rm -r pd_res_*")


#import sys
sys.stdout.flush()
print ('Content-Type: text/html\r\n')
#print('/r/n')
#os.system("rm -r pd_res_*")
if __name__ == '__main__':
	#print ('Content-Type: text/html\r\n')
	#print('/r/n')
	#print(pdb_id)
	#os.system("chmod +x pdpredict.py")
	print('<div style="height:400px">')
	print('</div>')
	print('<div style="vertical-align: middle;">')
	print("<center><div style='width:1000px;height:200px;border:solid black; font-size:20px;size = '+4''>")
	print("Feature calculations...")
	print("<br>")
	print("<center><div style='width:1000px;height:50px;border-top:solid black;font-size:15px;size = '+2''>")
	print('Please wait until the calculations are done')
	print('<br>')
	print(pdb_id,mutant)
	if pdb_id!="" and mutant!=""and pdb_file==b"" and mutfile==b"":
		pdb_id_up=pdb_id.upper()
		print("ppppppppppppppppppppppppppppppppppppppp")
		method=1
	elif pdb_id==""and pdb_file=="":
		print ('Content-Type: text/html\r\n')
		print('/r/n')
		print ('<html><head><title>Input error</title><body>PDB ID is missing</body></html>')
		exit()
	elif mutant=="" and mutfile=="":
		print ('Content-Type: text/html\r\n')
		print('/r/n')
		print ('<html><head><title>Input error</title><body>Provide the mutation</body></html>')
		exit()
	elif pdb_id=="" and pdb_file!=b"":
		fileitem=pdb_file
		method=2
		if fileitem.filename:
			fn = os.path.basename(fileitem.filename)
			open('tmp/' + fn, 'wb').write(fileitem.file.read())
	elif mutant=="" and mutfile!="":
		fileitem=mutfile
		method=3
		if fileitem.filename:
			fn = os.path.basename(fileitem.filename)
			open('tmp/' + fn, 'wb').write(fileitem.file.read())	
	else:
		print ('<html><head><title>Input error</title><body>No input found</body></html>')
		exit()
	
	pdb_id_up=pdb_id
	if method==1:
		dir_path = os.path.dirname(os.path.realpath(__file__))
		randname="pd_res_"+uuid.uuid4().hex
		path = os.path.join(dir_path, randname)
		os.mkdir(path,0o777)
		os.system("chmod -R 777 {}".format(path))
		#os.system("chmod 777 bin")
		#os.system("chmod +x bin")
		os.chdir(path)
		#print(os.listdir())
		#print(path)
		#shutil.copyfile("../1aay.pdb", "1aay.pdb")
		os.system("wget 'https://files.rcsb.org/download/{}.pdb1'".format(pdb_id_up))
		os.system("mv  "+pdb_id_up+".pdb1 "+pdb_id_up+".pdb")
		mut_list=mutant.split()
		#print(mut_list)
	elif method==2:
		dir_path = os.path.dirname(os.path.realpath(__file__))
		randname=uuid.uuid4().hex
		path = os.path.join(dir_path, randname)
		os.mkdir(path,0o777)
		os.system("chmod -R 777 {}".format(path))
		os.chdir(path)
		os.system("mv ../tmp/{} input.pdb".format(fn))
		pdb_id_up='input'
	elif method==3:
		dir_path = os.path.dirname(os.path.realpath(__file__))
		randname=uuid.uuid4().hex
		path = os.path.join(dir_path, randname)
		os.mkdir(path,0o777)
		os.system("chmod -R 777 {}".format(path))
		os.chdir(path)
		os.system("mv ../tmp/{} mutant.txt".format(fn))
		pdb_id_up='input'
		with open(r"mutant.txt") as file:
			mut_list=file.read().splitlines() 
		print(mut_list)
	
	print("</div></center></div>")
	print("<html>")
	print("<div style='width:100px;height:1px;overflow:hidden;'>")
	
	#shutil.copyfile("../clean_pdb.py", "clean_pdb.py")
	#os.system(r"python3 clean_pdb.py {}.pdb".format(pdb_id_up))
	#redirectURL = "%s/result.py" % randname
	#print ('Content-Type: text/html\r\n')
	#print ('<html>')
	#print ('  <head>')
	#print ('	<meta http-equiv="refresh" content="0;url=%s" />' % redirectURL)
	#print('	<title>You are going to be redirected</title>')
	class ProtSelect(Select):
		warnings.simplefilter('ignore', PDBConstructionWarning)
		warnings.simplefilter('ignore', FutureWarning)
		def accept_residue(self, residue):
			if not is_aa(residue, standard=True):
				res = residue.id[0]
				if not res == "W":
					return True
			else:
				return False
	class ProtSelect1(Select):
		def accept_residue(self, residue):
			warnings.simplefilter('ignore', PDBConstructionWarning)
			warnings.simplefilter('ignore', FutureWarning)
			if is_aa(residue, standard=True):
				return True
			else:
				return False
	
	parser = PDBParser()
	try:
		structure = parser.get_structure(pdb_id_up, pdb_id_up+".pdb")
	except:
		print("</div>")
		print("<center><div style='width:1000px;height:1000px;border-top:solid black;font-size:30px;color:red;size = '+2''>")
		print("Check if the PDB ID is correct")
		print("</div>")
		exit()
	from Bio.PDB import PDBParser, PDBIO, Select

	

	modell = structure[0]
	io = bpdb.PDBIO()
	io.set_structure(modell)
	modd=0
	io.save('dna_'+pdb_id_up+'.pdb', ProtSelect())
	io.save('apo_'+pdb_id_up+'.pdb', ProtSelect1())
	with open('dna_'+pdb_id_up+'.pdb') as filedna:
		with open('apo_'+pdb_id_up+'.pdb') as protfile:
			if len(filedna.readlines())<5 or len(protfile.readlines())<5:
				for modd in range(1,len(structure)):
					modell = structure[modd]
					io = bpdb.PDBIO()
					io.set_structure(modell)
					io.save('dna_'+pdb_id_up+'.pdb', ProtSelect())
					io.save('apo_'+pdb_id_up+'.pdb', ProtSelect1())
					if not len(filedna.readlines())<5 or not len(protfile.readlines())<5:
						break
	
	#shutil.copyfile("../foldx", "foldx")
	#shutil.copyfile("../rotabase.txt", "rotabase.txt")
	#shutil.copyfile("../naccess", "naccess")
	#shutil.copyfile("../bin","bin")
	shutil.copyfile("../style4.css", "style4.css")
	#shutil.copyfile("../clean_pdb.py", "clean_pdb.py")
	#os.system(r"python3 clean_pdb.py {}.pdb".format(pdb_id_up))
	shutil.copyfile("../index.txt", "index.txt")
	shutil.copyfile("../1asy_T124A.pdb", "1asy_T124A.pdb")
	shutil.copyfile("../footer.txt", "footer.txt")
	#shutil.copyfile("../dssp","dssp")
	shutil.copyfile("../rna_4vdrch.csv", "rna_4vdrch.csv")
	shutil.copyfile("../aa_20vdrch.csv", "aa_20vdrch.csv")
	shutil.copyfile("../potential_res.csv","potential_res.csv")
	#os.system("chmod -R +x {}".format("foldx"))
	#os.system("chmod -R +x {}".format("naccess"))
	#os.system("chmod -R 777 {}".format("naccess"))
	os.system("chmod -R +x {}".format("bin"))
	#os.system("chmod -R +x {}".format("dssp"))
	os.system("chmod -R 777 {}".format(path))
	d={"A":"ALA", "R":"ARG", "N":"ASN", "D":"ASP",
              "C":"CYS", "E":"GLU", "Q":"GLN", "G":"GLY",
              "H":"HIS", "I":"ILE", "L":"LEU", "K":"LYS",
              "M":"MET", "F":"PHE", "P":"PRO", "S":"SER",
              "T":"THR", "W":"TRP", "Y":"TYR", "V":"VAL"}
        
	print ('Content-Type: text/html\r\n')
	print('/r/n')
	d1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

	with open("result.txt","w") as resultout:
		for mut_ind in mut_list:
			print("</div>")
			print("<center><div style='width:1000px;height:50px;font-size:18px;color:black;size = '+2''>")
			print("Calculating the features for:"+mut_ind)
			print("</div>")
			print("<div style='width:100px;height:1px;overflow:hidden;'>")
			mut=mut_ind.split(":")
			mut=mut[1][0]+":"+mut[0]+":"+mut[1][1:-1]+":"+mut[1][-1]
			mut=mut.split(":")
			#print(mut)
			mut1=mut[0]+mut[2]+mut[3]
			chain=mut[1].split("-")
			#print(mut1,chain)
			parser = PDBParser()
			structure = parser.get_structure(pdb_id_up, pdb_id_up+".pdb")
			for chains in structure[modd]:
				if chains in chain:
					if len(chains)>=int(mut[2]):
						#print(list(chains))
						residue = chains[int(mut[2])]
						#print(residue.get_resname())
						#print(chains[(int(mut[2]))])
						#print(residue.get_resname())
						if d1[residue.get_resname()]!=mut[0]:
							print("</div>")
							print("<center><div style='width:1000px;height:1000px;border-top:solid black;font-size:30px;color:red;size = '+2''>")
							print("Wild-type position does not match")
							print("</div>")
							print("<div style='width:100px;height:1px;overflow:hidden;'>")
							exit()
			with open("pymoll_sample.py",'w') as file:
				file.write("import pymol\n")
				file.write("from pymol import cmd\n")
				file.write('cmd.wizard("mutagenesis")\n')
				file.write("cmd.load('"+pdb_id_up+".pdb')\n")
				file.write('cmd.get_wizard().set_mode("'+d[mut[3]]+'")\n')
				print(chain)
				file.write('for cha in '+str(chain)+':\n')
				file.write('	cha=cha.strip()\n')
				file.write('	cmd.get_wizard().do_select("chain "+cha+" and resid '+mut[2]+'")\n')
				file.write('	cmd.frame(1)\n')
				file.write('	cmd.get_wizard().apply()\n')
				file.write('cmd.save("'+path+'/'+pdb_id_up+"_"+mut1+'.pdb")\n')
			#os.chdir("/var/www/html/bioinfo2/cgi-bin/pdamutpred/Harini_Pymol/")
			print(path)
			print ('Content-Type: text/html\r\n')
			print('/r/n')
			#p2=subprocess.Popen("./foldx --command=AnalyseComplex --pdb-dir=../"+randname+" --pdb="+pdb_id_up+".pdb --complexWithDNA=True",  shell=True)
			p1=subprocess.Popen("/mnt/c/Users/ulian/Downloads/pra-mut-pred-git/Harini_Pymol/pymol/pymol -cq pymoll_sample.py",stdout=subprocess.PIPE, shell=True)
			p1.wait()
			class ProtSelect(Select):
				warnings.simplefilter('ignore', PDBConstructionWarning)
				warnings.simplefilter('ignore', FutureWarning)
				def accept_residue(self, residue):
					if not is_aa(residue, standard=True):
						res = residue.id[0]
						if not res == "W":
							return True
					else:
						return False
			class ProtSelect1(Select):
				def accept_residue(self, residue):
					warnings.simplefilter('ignore', PDBConstructionWarning)
					warnings.simplefilter('ignore', FutureWarning)
					if is_aa(residue, standard=True):
						return True
					else:
						return False
			os.system("chmod -R 777 "+pdb_id_up+"_"+mut1+".pdb")
			parser = PDBParser()
			#try:
			structure = parser.get_structure(pdb_id_up+"_"+mut1, pdb_id_up+"_"+mut1+".pdb")
			#except:
			#	print("<center><div style='width:1000px;height:1000px;border-top:solid black;font-size:30px;color:red;size = '+2''>")
			#	print("Check if the PDB ID is correct")
			#	print("</div>")
			#	exit()
			
			modell = structure[modd]
			io = bpdb.PDBIO()
			io.set_structure(modell)
			modd=0
			io.save('dna_'+pdb_id_up+"_"+mut1+".pdb", ProtSelect())
			io.save('apo_'+pdb_id_up+"_"+mut1+".pdb", ProtSelect1())
			#p1=subprocess.Popen("/var/www/html/bioinfo2/cgi-bin/pdamutpred/Harini_Pymol/PyMOL-3.0.0_appveyor664-Linux-x86_64-py310/pymol/pymol -cq pymoll_sample.py", shell=True)
			import sys
			sys.stdout.flush()
			print(os.getcwd())
			print(pdb_id_up)
#***********************************************************************************************************************************************************************************************
			os.chdir("../foldx5Linux64.tar_/")
			
			print(os.listdir)
			#os.system("chmod +x ./foldx_20231231")
			#sys.flush
			p2=subprocess.Popen("./foldx --command=AnalyseComplex --pdb-dir=../"+randname+" --pdb="+pdb_id_up+".pdb --output-dir=../"+randname+" --complexWithRNA=True", shell=True)
			p2.wait()
			fold_par=[]
			fold_val=[]
			fold_dict={}
			fold_dict={}
			try:
				for f in glob.glob('../'+randname+'/Interaction_'+pdb_id_up+'*.fxout'):
					#print(f)
					with open(f) as file:
						lis=file.readlines()
						#print(lis)
						par=lis[-2].split("\t")
						val=lis[-1].split("\t")
						for fo in range(len(val)):
								fold_dict[par[fo]]=val[fo]
					print(fold_dict)
				enegyw=float(fold_dict['Interaction Energy'])
			except:
				class NonHetSelect(Select):
					def accept_residue(self, residue):
						return 1 if residue.id[0] == " " else 0

				pdb = PDBParser().get_structure(pdb_id_up, "../"+randname+"/"+pdb_id_up+".pdb")
				io = PDBIO()
				io.set_structure(pdb)
				io.save("../"+randname+"/non"+pdb_id_up+".pdb", NonHetSelect())
				try:
					p2=subprocess.Popen("./foldx --command=AnalyseComplex --pdb-dir=../"+randname+" --pdb=non"+pdb_id_up+".pdb --output-dir=../"+randname+" --complexWithRNA=True",shell=True)
					p2.wait()
					for f in glob.glob('../'+randname+'/Interaction_'+pdb_id_up+'*.fxout'):
						with open(f) as file:
							lis=file.readlines()
							#print(lis)
							par=lis[-2].split("\t")
							val=lis[-1].split("\t")
							for fo in range(len(val)):
									fold_dict[par[fo]]=val[fo]
						print(fold_dict)
					enegyw=float(fold_dict['Interaction Energy'])
				except:
					print("check PDB ID, problem with foldx")
					exit()
			p2=subprocess.Popen("./foldx --command=AnalyseComplex --pdb-dir=../"+randname+" --pdb="+pdb_id_up+"_"+mut1+".pdb --output-dir=../"+randname+" --complexWithRNA=True",shell=True)
			p2.wait()
			fold_par=[]
			fold_val=[]
			fold_dict={}
			fold_dict={}
			try:
				for f in glob.glob('../'+randname+'/Interaction_'+pdb_id_up+"_"+mut1+'*.fxout'):
					#print(f)
					with open(f) as file:
						lis=file.readlines()
						#print(lis)
						par=lis[-2].split("\t")
						val=lis[-1].split("\t")
						for fo in range(len(val)):
								fold_dict[par[fo]]=val[fo]
					print(fold_dict)
				enegym=float(fold_dict['Interaction Energy'])
			except:
				class NonHetSelect(Select):
					def accept_residue(self, residue):
						return 1 if residue.id[0] == " " else 0

				pdb = PDBParser().get_structure(pdb_id_up, "../"+randname+"/"+pdb_id_up+"_"+mut1+".pdb")
				io = PDBIO()
				io.set_structure(pdb)
				io.save("../"+randname+"/non"+pdb_id_up+"_"+mut1+".pdb", NonHetSelect())
				try:
					p2=subprocess.Popen("./foldx --command=AnalyseComplex --pdb-dir=../"+randname+" --pdb=non"+pdb_id_up+"_"+mut1+".pdb --output-dir=../"+randname+" --complexWithRNA=True",shell=True)
					p2.wait()	
					for f in glob.glob('../'+randname+'/Interaction_'+pdb_id_up+"_"+mut1+'*.fxout'):
					#print(f)
						with open(f) as file:
							lis=file.readlines()
							#print(lis)
							par=lis[-2].split("\t")
							val=lis[-1].split("\t")
							for fo in range(len(val)):
									fold_dict[par[fo]]=val[fo]
						print(fold_dict)
					enegym=float(fold_dict['Interaction Energy'])
				except:
					print("check PDB ID, problem with foldx")
					exit()
			print(enegym,enegyw)
			os.chdir("../"+randname)
#***********************************************************************************************************************************************************************************************'''
			prot_chain=[]
			dna_chain=[]
			with open("apo_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							prot_chain.append(x)
			with open("dna_"+pdb_id_up+".pdb") as file:
				for rows in file.readlines():
					if rows[0:4]=='ATOM':
						x=rows[21].strip()
						if not x=='':
							dna_chain.append(x)	
			dna_chain=list(set(dna_chain))
			prot_chain=list(set(prot_chain))
#***********************************************************************************************************************************************************************************************'''
			os.system(r"chmod -x ../stride/stride")
			os.system(r"../stride/stride -h apo_"+pdb_id_up+"_"+mut1+".pdb -fmut_apo_"+pdb_id_up+"_"+mut1+"_stride.out")
			os.system(r"../stride/stride -h apo_"+pdb_id_up+".pdb -fapo_"+pdb_id_up+"_stride.out")
			rr=prot_chain
			print(prot_chain)
			mut_area=[]
			wild_area=[]
			with open(r"mut_apo_"+pdb_id_up+"_"+mut1+"_stride.out") as file:
				for j in file.readlines():
					#print(x[i])
					if j[0:3]=="ASG":
						if j[5:8] in d1:
							if j[11:15][-1].isalpha():
								y=j[11:14]
							else:
								y=j[11:15]
							#print(mut[2],mut[3])
							#print(j)
							#print(rr)
							#print(j)
							if d1[j[5:8]]==mut[3] and j[9] in rr and int(y)==int(mut[2]):
								mut_area.append(float(j[64:69]))
			with open(r"apo_"+pdb_id_up+"_stride.out") as file:
				for j in file.readlines():
					#print(x[i])
					if j[0:3]=="ASG":
						if j[5:8] in d1:
							if j[11:15][-1].isalpha():
								y=j[11:14]
							else:
								y=j[11:15]
							#print(mut[1:-1])
							#print(j)
							#print(rr)
							#print(j)
							if d1[j[5:8]]==mut[0] and j[9] in rr and int(y)==int(mut[2]):
								wild_area.append(float(j[64:69]))
			print(wild_area,mut_area)
			diff_area=sum(mut_area)/len(mut_area)-sum(wild_area)/len(wild_area)
			print(diff_area)
			#print(ii)	
#***********************************************************************************************************************************************************************************************'''
			def map_sift(filename):
				with open(filename) as file:
					pdb=filename.split(".")[0]
					#print(filename)
					data={}
					ins,un,sc,ins1,ss=1,1,1,0,0
					start_end,segid_number,residue,single,name,start_end,segid,pdb_pos,chain,uniprot,uni_name,scop=[],[],[],[],[],[],[],[],[],[],[],[]
					for i in file.readlines():
						if re.match('(.*)entityId="(.*)"',i):
							x=re.match('(.*)entityId="(.*)"',i).group(2)
							#chain=x
						if re.match('(.*)segment (.*) start="(.*)" end="(.*)"',i):
						    y=re.match('(.*)segment (.*) start="(.*)" end="(.*)"',i)
						    #print(y.group(3))
						    start_end_m=(str(y.group(3))+"-"+str(y.group(4)))
						    segid_number=y.group(2)
						if re.match('(.*)residue dbSource="(.*)" dbCoordSys="(.*)" dbResNum="(.*)" dbResName="(.*)"(.*)',i):
						    res=re.match('(.*)residue dbSource="(.*)" dbCoordSys="(.*)" dbResNum="(.*)" dbResName="(.*)"(.*)',i)
						    residue.append(res.group(4))
						    name.append(res.group(5))
						    start_end.append(start_end_m)
						    segid.append(segid_number)
						    ins=0
						if re.match('(.*)crossRefDb dbSource="PDB" dbCoordSys="(.*)" dbAccessionId="(.*)" dbResNum="(.*)" dbResName="(.*)" dbChainId="(.*)"/>',i):
						    pb=re.match('(.*)crossRefDb dbSource="PDB" dbCoordSys="(.*)" dbAccessionId="(.*)" dbResNum="(.*)" dbResName="(.*)" dbChainId="(.*)"/>',i)
						    pdb_pos.append(pb.group(4))
						    chain.append(pb.group(6))
						if re.match('(.*)crossRefDb dbSource="UniProt" dbCoordSys="(.*)" dbAccessionId="(.*)" dbResNum="(.*)" dbResName="(.*)"(.*)',i):
						    uni=re.match('(.*)crossRefDb dbSource="UniProt" dbCoordSys="(.*)" dbAccessionId="(.*)" dbResNum="(.*)" dbResName="(.*)"(.*)',i)
						    uniprot.append(uni.group(4))
						    uni_name.append(uni.group(3))
						    single.append(uni.group(5))
						    un=0
						if re.match("(.*)</residue>",i):
							if un==1:
								if ins==0:
									uniprot.append("-")
									uni_name.append("-")
									single.append("-")
							ins,un,sc,ins1,ss=1,1,1,0,0
						#print(len(start_end),len(segid),len(residue),len(name),len(start_end),len(segid),len(pdb_pos),len(chain),len(uniprot),len(uni_name),len(scop))
				data={'chain':chain,'seg':segid,'uniprot_id':uni_name,'start_end':start_end,'residue_no':residue,'name':name,'residue':single,'pdb_position':pdb_pos,'uniprot_pos':uniprot}
				df=pd.DataFrame(data)
				df.to_csv(r"sifts_mapped_mutation_"+pdb+".csv")
	
			os.system("wget https://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/"+pdb_id_up[1:-1]+"/"+pdb_id_up+".xml.gz")
			os.system("chmod -R 777 "+pdb_id_up+".xml.gz")
			os.system("gzip -d "+pdb_id_up+".xml.gz")
			map_sift(pdb_id_up+".xml")
			pop={}
			df2=pd.read_csv(r"sifts_mapped_mutation_"+pdb_id_up+".csv")
			for k in mut[1].split("-"):
				k=k.strip()
				df1=df2.loc[df2["chain"]==k]
				df1=df1.dropna()
				#print(df1)
				df1["new_residue"]=range(1,len(df1)+1)
				pos=mut[2]
				wild=mut[0]
				#print(mut)
				df1=df1.loc[(df1["pdb_position"]==int(pos))&(df1["residue"]==wild)]
				pop[k]=df1["new_residue"].values
			parser = PDBParser(QUIET=True)
			structure = parser.get_structure('struct', pdb_id_up+".pdb")    

			# iterate each model, chain, and residue
			# printing out the sequence for each chain
			#cha=[p.strip() for p in k.split(":")[1].split(",")]
			#print(cha)
			with open(pdb_id_up+"_list.txt","a+") as file1:
			#for model in structure[0]:
			     for chain1 in structure[0]:
			     	#print(chain)
			     	if chain1.id in chain:
				     	seq = []
				     	for residue in chain1:
				     		if residue.resname in d1:
					     		seq.append(d1[residue.resname])
			     		print('>'+pdb_id_up+'\n',''.join(seq),file=open(pdb_id_up+"_"+chain1.id+".fasta","a"))
			     		file1.write(pdb_id_up+"_"+chain1.id+".fasta")
			     		file1.write("\n")
			f=open(pdb_id_up+"_list.txt").readlines()
			import os
			for i in f:
				i=i.rstrip().replace(".fasta","")
				#print (i)
				#os.chdir("/mnt/c/Users/ulian/Downloads/pra-mut-pred-git/ncbi-blast-2.16.0+/bin")
				os.system("/mnt/c/Users/ulian/Downloads/pra-mut-pred-git/ncbi-blast-2.16.0+/bin/psiblast -query {0}.fasta -db ../uniprot_db/uniprot_sprot.fasta -out {0}.out -num_iterations 3 -num_threads 8 -out_ascii_pssm {0}_out.pssm".format(i))
				os.system("/mnt/c/Users/ulian/Downloads/pra-mut-pred-git/ncbi-blast-2.16.0+/bin/blastp -query {0}.fasta -db ../uniprot_db/uniprot_sprot.fasta -num_threads 8 -outfmt 4 -out {0}.out -max_target_seqs 100".format(i))
			data=['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
			f=open(pdb_id_up+"_list.txt").readlines()
			#mutlist=open("muts_list.txt").readlines()
			info_per=[]
			print ("   pssm_mut pssm_mut_weighted pssm_diff pssm_weight_diff info_per_pos rel_weight_to_pseudo")
			for ch in mut[1].split("-"):
				ch=ch.strip()
				#pdb=df["pdb_selected_new"][i].split("_")[0]
				print(mut)
				ip=mut[0]+ch+mut[2]+mut[3]
				wt=mut[0]
				mutation=mut[-1]
				number=pop[ch]
				g=open(pdb_id_up+"_"+ch+"_out.pssm").readlines()[3:-6]
				for j in g:
					j=j.rstrip().split()
					if int(j[0])==number:
						info_per.append(float(j[-2]))
						print (ch,ip[0],j[data.index(mutation)+2],j[data.index(mutation)+22],float(j[data.index(mutation)+2])-float(j[data.index(wt)+2]),float(j[data.index(mutation)+22])-float(j[data.index(wt)+22]),j[-2],j[-1])
			#print(oo)
			infor_per_pos=sum(info_per)/len(info_per)
			
#***********************************************************************************************************************************************************************************************'''

			#def calculate_half_sphere_exposure(structure, radius=10.0):
			    # Initialize a dictionary to store half-sphere exposure for each residue
			hse_values = {}
			radius=10.0
			#pdb_file = "1fos.pdb"
			parser = PDB.PDBParser()
			structure = parser.get_structure("protein",pdb_id_up+".pdb")
			hse=[]
			#print(x[pp])
			for chains in structure[0]:
				if chains.id in chain:
					print(chains)
					residues = list(chains)
					total_residues = len(residues)
					#res=prot_res.resname
					print(mut)
					for i,residue_i in enumerate(chains):
						if (int(residue_i.get_full_id()[3][1])==int(mut[2]) and d1[residue_i.resname]==mut[0]):
							print(residue_i)
							atom_list = list(residue_i.get_atoms())
							atom_list=[k.id for k in atom_list]
							if "CA" in atom_list:
								ca_coords_i = residue_i['CA'].get_coord()
								exposed_ca_atoms = 0
								for j, residue_j in enumerate(chains):
									if i != j:
										try:
											ca_coords_j = residue_j['CA'].get_coord()
											distance = np.linalg.norm(ca_coords_i - ca_coords_j)
											if distance <= radius:
												exposed_ca_atoms += 1
										except KeyError:
											continue
								hse_value = exposed_ca_atoms / total_residues
								#node_name = residue_i
								hse.append(hse_value)
			hse_w = sum(hse)/len(hse)
			print(hse_w)
			os.system("../x3dna-dssr --get-hbond --json -i="+pdb_id_up+".pdb -o=hbond_"+pdb_id_up+".out")
			df1=pd.read_json(r"hbond_"+pdb_id_up+".out")
			df1=pd.DataFrame(dict(df1["hbonds"]))
			df1=df1.T
			df1.to_csv(r"csv_hbond_"+pdb_id_up+".out")
			df=pd.read_csv(r"csv_hbond_"+pdb_id_up+".out")
			nt_aa=len(df.loc[df["residue_pair"]=="nt:aa"])
			print(nt_aa)
#***********************************************************************************************************************************************************************************************'''
			cv_res=[]
			radius=6.0
			structure = parser.get_structure("protein",pdb_id_up+".pdb")
			for chains in structure[0]:
				if chains.id in chain:
					for i,residue_i in enumerate(chains):
						cv=[]
						if (int(residue_i.get_full_id()[3][1])==int(mut[2]) and d1[residue_i.resname]==mut[0]):
							for atoms in residue_i:
								vector=[]
								ca_coords_i = atoms.get_coord()
								for chains in structure[0]:
									#if chains.id in df["prot_chain1"][pp]:
									for j, residue_j in enumerate(chains):
										if i != j:
											for atom1 in residue_j:
												ca_coords_j = atom1.get_coord()
												distance = np.linalg.norm(ca_coords_i - ca_coords_j)
												if distance <= radius:
													#print(atom1.get_coord(), atoms.get_coord())
													#print(atom1.get_coord()-atoms.get_coord())
													distance = np.linalg.norm(atom1.get_coord() - atoms.get_coord())
													vector.append((atom1.get_coord()-atoms.get_coord())/distance)
								ni=len(vector)
								sumation=np.sum(np.array(vector),axis=0)
								magnitude=np.sqrt(np.sum(np.square(sumation)))
								#print(ni,magnitude)
								if ni!=0:
									cv.append(1-(magnitude/ni))
								#print(cv)
								#cv = [x for x in cv if ~np.isnan(x)]
							cv_res.append(sum(cv)/len(cv))
							#print(cv_res)
			print(len(cv_res))
			cv_res_w=sum(cv_res)/len(cv_res)
			structure = parser.get_structure("protein",pdb_id_up+"_"+mut1+".pdb")
			cv_res=[]
			for chains in structure[0]:
				if chains.id in chain:
					for i,residue_i in enumerate(chains):
						cv=[]
						if (int(residue_i.get_full_id()[3][1])==int(mut[2]) and d1[residue_i.resname]==mut[-1]):
							for atoms in residue_i:
								vector=[]
								ca_coords_i = atoms.get_coord()
								for chains in structure[0]:
									#if chains.id in df["prot_chain1"][pp]:
									for j, residue_j in enumerate(chains):
										if i != j:
											for atom1 in residue_j:
												ca_coords_j = atom1.get_coord()
												distance = np.linalg.norm(ca_coords_i - ca_coords_j)
												if distance <= radius:
													#print(atom1.get_coord(), atoms.get_coord())
													#print(atom1.get_coord()-atoms.get_coord())
													distance = np.linalg.norm(atom1.get_coord() - atoms.get_coord())
													vector.append((atom1.get_coord()-atoms.get_coord())/distance)
								ni=len(vector)
								sumation=np.sum(np.array(vector),axis=0)
								magnitude=np.sqrt(np.sum(np.square(sumation)))
								#print(ni,magnitude)
								if ni!=0:
									cv.append(1-(magnitude/ni))
								#print(cv)
								#cv = [x for x in cv if ~np.isnan(x)]
							cv_res.append(sum(cv)/len(cv))
							#print(cv_res)
			print(len(cv_res))
			cv_res_m=sum(cv_res)/len(cv_res)
			print(cv_res_m,cv_res_w)
			cv_diff=cv_res_m-cv_res_w
			print(cv_diff)
			import urllib.request 
			class Protein_RNA_ineractions:
				def __init__(self, pdb_file,prot_chain='A',rna_chain='B'):
					self.pdb_file = pdb_file
					self.prot_chain = prot_chain
					self.rna_chain = rna_chain
					self.pattern ='^ATOM.{16}'
				def f6_proteinRNAcontact2(self,pdb,ppp,rrr,dist_cutoff=3.5):
					warnings.simplefilter('ignore', PDBConstructionWarning)
					warnings.simplefilter('ignore', FutureWarning)
					int_df1 = pd.DataFrame()
					flag = 0
					parser = PDB.PDBParser()
					structure = parser.get_structure("pdb", pdb)
					model = structure[modd]
					prot_chain = model[self.prot_chain] 
					rna_chain = model[self.rna_chain]
					nal1 = ['A','G','C','T','U']
					pal1 = ['ALA','ARG','ASN','ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
					c1 = [str(x).split()[1] for x in list(prot_chain.get_residues())]
					c2 = [str(x).split()[1] for x in list(rna_chain.get_residues())]
					checkp = set(c1).intersection(set(nal1))
					checkn = set(c2).intersection(set(pal1))
					if (len(checkp) > 0):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See protein set "+str(checkp))
						flag = 1
						#return int_df1, flag
					if (len(checkn) >= 2):
						print("check the chains of pdb file "+self.pdb_file+" error in protein or RNA molecule\n See na set "+(str(checkn)))
						flag = 1
						#return int_df1, flag
					v=1
					##print("")
					for prot_res in prot_chain:
						if prot_res.resname in d1:
							for prot_atoms in prot_res:
								for rna_res in rna_chain:
									rna_resname = rna_res.resname
									for rna_atoms in rna_res:
										distance = prot_atoms-rna_atoms
										##print("")
										if (distance<= 6):
											dict1 = {'distance':distance, 'na_atm':rna_atoms.get_full_id()[4][0], 'na_atmno':rna_atoms.get_serial_number(),'na_res':rna_res.resname, 'na_resno':rna_atoms.get_full_id()[3][1], 'na_coord':rna_atoms.get_coord(), 'prot_atm':prot_atoms.get_full_id()[4][0],'prot_atmno':prot_atoms.get_serial_number(), 'prot_res':prot_res.resname, 'prot_resno':prot_atoms.get_full_id()[3][1],'prot_coord':prot_atoms.get_coord()}
											int_df1 = int_df1.append(dict1,ignore_index=True)
											v=0
											
						if v==0:
							df_inter=pd.DataFrame(int_df1)
						#b1 = [x.strip() in nal1 for x in int_df1['na_res']]
						#temp1 = int_df1[b1]
						#b2 = [x.strip() in pal1 for x in temp1['prot_res']]
						#df_inter_n = temp1[b2]
						
					else:
						df_inter=pd.DataFrame(int_df1)
					#print(df_inter)
					return df_inter
				def f8_interaction_type (self, df_inter):
					int_dict = {'CO':0,'OC':0,'NO':0,'ON':0}
					ppppp=['CO','OC','NO','ON','CP','NP','NN','CN','CC','NC','OO','OP','SP','SO','SC','SN']
					print(df_inter)
					temp_ptr = open('temp.txt','w')
					for index, dfrow in df_inter.iterrows():
						temp_ptr.write('\n'+dfrow["prot_atm"]+'\t'+dfrow["na_atm"])
						in_t = dfrow["prot_atm"][0:1]+dfrow["na_atm"][0:1]
						if in_t in int_dict:
							int_dict[in_t] += 1
						else:
							if not 'H' in in_t and in_t in ppppp:
								int_dict.update({in_t:1})
					temp_ptr.close()
					return int_dict
			
			aa_param = pd.read_csv('aa_20vdrch.csv')
			na_param = pd.read_csv('rna_4vdrch.csv')
			res_bind_count=[]	
			atom_count=[]
			scsc=[]
			catom=0
			bind=[]
			tot=[]
			op_count,sc_count,nn_count,on_count,co_count=[],[],[],[],[]
			scm_count,nom_count,onm_count=[],[],[]
			totm=[]
			for i in prot_chain:
				for j in dna_chain:
					inst2 = Protein_RNA_ineractions(pdb_id_up+'.pdb', prot_chain=i,rna_chain=j)
					df_inter=inst2.f6_proteinRNAcontact2(pdb=pdb_id_up+'.pdb',ppp=i,rrr=j)
					df_inter = pd.DataFrame(df_inter)
					##print("")
					if len(df_inter.columns)>3:
						df_inter.to_csv(pdb_id_up+"_"+i[0]+"_"+j[0]+"_atom_interaction.csv")
						###print("")
						atoms_involve = inst2.f8_interaction_type(df_inter)
						###print("")
						
						df_at=pd.Series(atoms_involve).to_frame()
						df_at.columns=['count']
						df_at.to_csv(pdb_id_up+"_"+i+"_"+j+"_atom_count.csv")
						tot.append(df_at['count'].sum())
						#OP', 'NC1', 'Incl', 'CN1', 'OC1']
						if 'OP' in df_at.index:
							op_count.append(int(df_at.loc[['OP']].values))#/df_at['count'].sum())	
						if 'SC' in df_at.index:
							sc_count.append(int(df_at.loc[['SC']].values))#/df_at['count'].sum())
						if 'NN' in df_at.index:
							nn_count.append(int(df_at.loc[['NN']].values))#/df_at['count'].sum())
						if 'ON' in df_at.index:
							on_count.append(int(df_at.loc[['ON']].values))#/df_at['count'].sum())
						if 'CO' in df_at.index:
							co_count.append(int(df_at.loc[['CO']].values))#/df_at['count'].sum())
						##print("")
					df_inter=inst2.f6_proteinRNAcontact2(pdb=pdb_id_up+"_"+mut1+'.pdb',ppp=i,rrr=j)
					df_inter = pd.DataFrame(df_inter)
					##print("")
					if len(df_inter.columns)>3:
						df_inter.to_csv("mut_"+pdb_id_up+"_"+mut1+"_"+i[0]+"_"+j[0]+"_atom_interaction.csv")
						##print("")
						atoms_involve = inst2.f8_interaction_type(df_inter)
						##print("")
						df_at=pd.Series(atoms_involve).to_frame()
						df_at.columns=['count']
						df_at.to_csv("mut_"+pdb_id_up+"_"+mut1+"_"+i+"_"+j+"_atom_count.csv")
						totm.append(df_at['count'].sum())
						if 'SC' in df_at.index:
							scm_count.append(int(df_at.loc[['SO']].values))#/df_at['count'].sum())
						if 'NO' in df_at.index:
							nom_count.append(int(df_at.loc[['NO']].values))#/df_at['count'].sum())
						if 'ON' in df_at.index:
							onm_count.append(int(df_at.loc[['ON']].values))#/df_at['count'].sum())
						##print("")
			print(sum(op_count),sc_count,nn_count,on_count,co_count,scm_count,nom_count,onm_count)
			print(sum(op_count)/sum(tot),sum(sc_count)/sum(tot),sum(nn_count),sum(on_count)/sum(tot),sum(co_count)/sum(tot),sum(scm_count)/sum(totm),sum(nom_count)/sum(totm),sum(onm_count)/sum(totm))	
			#print(ii)
			print(sum(totm))
			opw=sum(op_count)/sum(tot)
			#print(sum(scm_count)/sum(totm))
			#print(sum(sc_count)/sum(tot))
			scmut=sum(scm_count)/sum(totm)
			scw=sum(sc_count)/sum(tot)
			delsc=scmut-scw
			#print(scmut)
			#print(scw)
			nn=sum(nn_count)
			inter_diff=enegym-enegyw
			nom=sum(nom_count)/sum(totm)
			on=(sum(onm_count)/sum(totm))-(sum(on_count)/sum(tot))
			cow=sum(co_count)/sum(tot)
			#print(opw,delsc,nn,inter_diff,nom,on,cow)
			#print(sum(sc_count),sum(scm_count),tot,totm)
			#print(ii)
			from collections import Counter
#***********************************************************************************************************************************************************************************************'''
			df1=pd.read_csv(r"../res_potential_final.csv")
			df1=df1.set_index('aa')
			prefix=pdb_id_up
			pdb=[]
			j=0
			for file in os.listdir(): 
				if file.startswith(prefix) and "_atom_interaction.csv" in file and ".csv" in file and not "mut" in file: 
					print(prefix)
					print(file)
					#natype=df.loc[df["pdb"]==prefix]["na"].values[0]
					list_of_dfs= [pd.concat( [pd.read_csv (file)],ignore_index=True)]
					#print(list_of_dfs)
					list_of_dfs = [df for df in list_of_dfs if len(df.columns)>5]
					#print(list_of_dfs)
					if len(list_of_dfs)>0:
						if not prefix in pdb:
							pdb.append(prefix)
							final=pd.concat(list_of_dfs)
							#print(final)
							final.drop(final.loc[final['prot_res']=="HOH"].index, inplace=True)
							final.drop(final.loc[final['na_res']=="HOH"].index, inplace=True)
							#final['prot_atm1'] = final['prot_atm'].astype(str).str[0]
							#final['na_atm1'] = final['na_atm'].astype(str).str[0]
							if "prot_res" in final.columns:
								final1=final[['prot_res','prot_resno','na_res','na_resno']]
								final1.drop_duplicates()
								print(final1)
								final1['combined'] = final1[['prot_res','na_res']].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
								print(final1)
								c = dict(Counter(final1['combined'].tolist()))
								j=0
								print(c)
								for i in c:
									k=i.split("_")
									if k[1] in df1.columns:
										j=j+(c[i]*float(df1.loc[k[0]][k[1]]))
										print(j)
									
								res_pot=j	
								print(j)
				
			print(res_pot)	
			#print(res_pot_mut)
			from sklearn.pipeline import make_pipeline
			from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
			from scipy.stats import pearsonr, spearmanr
			import joblib
			from sklearn.svm import SVR
			from sklearn.preprocessing import StandardScaler
			dataset_path = '../new_rna_train5.csv'

			dataframex = pd.read_csv(dataset_path)
			dataframex=dataframex.loc[dataframex["pdb_selected_new"]!="1urn_E19F"]
			dataframex=dataframex.loc[dataframex["pdb_selected_new"]!="1urn_E19S"]
			dataframex=dataframex.loc[dataframex["pdb_selected_new"]!="1urn_F19E"]
			dataframex=dataframex.loc[dataframex["pdb_selected_new"]!="1urn_S19E"]

			combo=['CO_wild_1', 'Interaction Energy_diff', 'NN', 'NO_mut_1', 'ON_wild_1_diff', 'OP_wild_1', 'SC_wild_1_diff', 'cv_diff', 'diff_area', 'half_sphere_exposure_wild', 'info_per_pos', 'nt:aa_x', 'res_pot']
			# Assuming 'X' contains your features, 'exp' is the target variable
			X = dataframex[list(combo)].values
			y = dataframex.iloc[:, 29].values
		
			# Create and train the GradientBoostingRegressor model with tuning conditions and a random seed
			model = make_pipeline(StandardScaler(), SVR(C=5, degree= 8, epsilon=0.2, gamma= 'scale',kernel='rbf'))
			model.fit(X, y)

			# Make predictions on the same dataset
			predictions = model.predict(X)
			pdb=dataframex["pdb"].tolist()
			# Calculate Mean Absolute Error (MAE)
			mae = mean_absolute_error(y, predictions)
			print(f'Mean Absolute Error (MAE): {mae}')

			# Calculate Pearson correlation and its p-value
			pearson_corr, p_value_pearson = pearsonr(y, predictions)
			#for kk in range(len(y)):
			#	print(y[kk], predictions[kk],pdb[kk])
			#print(ii)
			df=pd.DataFrame({'pdb':dataframex["pdb_selected_new"].tolist(),'y':y,'yp':predictions})
			df.to_csv(r"prediction.csv")
			print(f'Pearson Correlation: {round(pearson_corr, 3)}')
			print(f'P-value (Pearson): {p_value_pearson}')

			# Save the trained model to a file
			joblib.dump(model, 'trained_model.pkl')

			# Load the trained model from the file
			loaded_model = joblib.load('trained_model.pkl')

			# Function to make predictions
			def predict(features):
				
				return loaded_model.predict([features])
			
			# Example usage
			print(cow,inter_diff,nn,nom,on,opw,delsc,cv_diff,diff_area,hse_w,infor_per_pos,nt_aa,res_pot)
			#print(ii)
			new_features = [cow,inter_diff,nn,nom,on,opw,delsc,cv_diff,diff_area,hse_w,infor_per_pos,nt_aa,res_pot]#       delsc, nn, inter_diff, infor_per_pos, diff_area, nom, on, hse_w, nt_aa, cv_diff, cow,res_pot,res_pot_mut]
			#new_features = [opw, delsc, nn, inter_diff, infor_per_pos, diff_area, nom, on, hse_w, nt_aa, cv_diff, cow,res_pot,res_pot_mut]
			predictions = predict(new_features)
			predictions = predictions[0]
			predval = f"{predictions:.4f}"
			
			# Determine the nature of affinity change
			if float(predval) >= 0:
				affinity_change = "Decrease"
			else:
				affinity_change = "Increase"	

			#predictions="%.2f" % predictions
			resultout.write(mut_ind)
			resultout.write("\t")
			resultout.write(str(predval))
			resultout.write("\t")
			if float(predval)>0:
				resultout.write("Decreasing affinity")
			if float(predval)<0:
				resultout.write("Increasing affinity")
			if float(predval)==0:
				resultout.write("No change in affinity")	
			resultout.write(str("\n"))
			for i in resultout.readlines():
				print("Mutation:"+i.split("\t")[0]+", Predicted ∆∆G (kcal/mol):"+i.split("\t")[1]+", Effect of mutation:"+i.split("\t")[2])
print("</div>")	

'''
timetaken=timeit.default_timer() - start
print("check")
os.system("chmod +x result.py")
with open("result.py", "w") as polyout:
	print("ii")
	polyout.write("""#!/opt/websites/anaconda/envs/tf37/bin/python\nimport cgi\nimport cgitb; cgitb.enable()\nprint ('Content-Type: text/html\\r\\n')\n""")
	g=open("index.txt").readlines()
	for gg in g:
		gg=gg.rstrip()
		polyout.write("""print (\"\"\"{}\"\"\")""".format(gg))
		polyout.write("\n")
	polyout.write("""print ('<font face="Times New Roman" ><table id="customers">')""")
	polyout.write("\n")
	polyout.write("print ('<tr>')")
	polyout.write("\n")
	polyout.write("print ('<td><b>Time taken</b></td><td>{:.2f} seconds</td></td>')".format(timetaken))
	polyout.write("\n")
	polyout.write("print ('</tr>')")
	polyout.write("\n")
	polyout.write("print ('<tr>')")
	polyout.write("\n")
	polyout.write("print ('<td><b>PDB ID</b></td><td>{}</td>')".format(pdb_id_up))
	polyout.write("\n")
	polyout.write("print ('</tr>')")
	polyout.write("\n")
	polyout.write("print ('</table></font>')")
	polyout.write("\n")
	polyout.write("print ('<br/>')")
	polyout.write("\n")
	polyout.write("print ('<br/>')")
	polyout.write("\n")
	polyout.write("""print ('<font face="Times New Roman" ><table id="customers">')""")
	polyout.write("\n")
	polyout.write("print ('<tr>')")
	polyout.write("\n")
	polyout.write("print ('<td><b>Mutation</b></td><td><b>Predicted &Delta;&Delta;G (kcal/mol)</b></td><td><b>Effect of mutation</b></td>')")
	polyout.write("\n")
	polyout.write("print ('</tr>')")
	polyout.write("\n")
	f=open("result.txt").readlines()
	for mut in range(len(f)):
		mut=f[mut].split()
		polyout.write("print ('<tr>')")
		polyout.write("\n")
		polyout.write("print ('<td>{}</td><td>{}</td><td>{}</td>')".format(mut[0].rstrip(),mut[1].rstrip(),mut[2].rstrip()))
		polyout.write("\n")
		polyout.write("print ('</tr>')")
		polyout.write("\n")
	
	polyout.write("print ('</table></font>')")
	polyout.write("\n")
	h=open("footer.txt").readlines()
	for hh in h:
		hh=hh.rstrip()
		polyout.write("""print (\"\"\"{}\"\"\")""".format(hh))
		polyout.write("\n")
		print("check")
print("</div>")
'''
os.system("chmod +x result.py")
os.system("ls > remfile")

'''
redirectURL = "/cgi-bin/%s/result.py" % randname
print ('Content-Type: text/html\r\n')
print ('<html>')
print ('  <head>')
print ('	<meta http-equiv="refresh" content="0;url=%s" />' % redirectURL)
print ('	<title>You are going to be redirected</title>')

'''

redirectURL = "%s/result.py" % randname
print ('Content-Type: text/html\r\n')
print ('<html>')
print ('  <head>')
print ('	<meta http-equiv="refresh" content="0;url=%s" />' % redirectURL)
print('	<title>You are going to be redirected</title>')

remfile=open("remfile").readlines()
for rf in remfile:
	rf=rf.rstrip()
	if rf =='molecules':
		os.rmdir('{}'.format(rf))
	elif rf =='autodock':
		shutil.rmtree('autodock')
	elif rf not in ['result.txt','result.py',"style4.css"]:
		os.remove('{}'.format(rf))

#print ('	Redirecting... <a href="%s">Click here if you are not redirected</a>' % redirectURL)
print ('</body>')
print ('</html>')

########## common for all the classifications - end
