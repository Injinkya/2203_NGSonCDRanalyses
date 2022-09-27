#!/usr/bin/python3

import os
import sys


#Variables here
#################################
##YOU NEED TO EDIT DIRPATH#######
#################################
dirpath = sys.argv[1]

#Input Parental H3 and Parental L3 sequence
ParH3 = sys.argv[2]
ParL3 = sys.argv[3]
if ParL3 == "NONE":
	ParL3=""

#Reads with frequency <= nThr will be ignored in the output
nThr = int(sys.argv[4])

#opens files to read and write
for fName in os.listdir(dirpath):
	if (fName.endswith(".fa")) and ("_no.fa" not in fName):
		fFasta = open(os.path.join(dirpath,fName),'r')
		fProtHtml = open(os.path.join(dirpath,fName[:-3]+"_"+str(nThr)+".html"),'w')
		fProtTxt = open(os.path.join(dirpath,fName[:-3]+"_"+str(nThr)+".txt"),'w')
		fProtHeat = open(os.path.join(dirpath,fName[:-3]+"_"+str(nThr)+"_heat.txt"),'w')

#make freq dictionary dCount[sequence]=number of reads
def freq_dict(dCount,sKey):
        try:
                dCount[sKey] += 1
        except KeyError:
                dCount[sKey] =1
        return dCount

#protein translate
def translate(seq):
	
	table = {
		'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
		'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
		'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',				
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
		'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
		'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
		'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
		'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
		'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
		'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
		'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
		'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
		'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
		'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
	}
	protein =""
	if len(seq)%3 == 0:
		for i in range(0, len(seq), 3):
			codon = seq[i:i + 3]
			protein+= table[codon]
	else:
		pass

	return protein

#count mutations in each position
def count_mutation(seq,Par1,Par2):
	smut=""
	for i in range(0,len(Par1)):
		if (seq[i]==Par1[i]):
			pass
		else:
			smut=smut+str(i+1)+seq[i]+" "
	if len(Par2)!=0:
		for i in range(len(Par1)+1,len(seq)):
			if (seq[i]==Par2[i-1-len(Par1)]):
				pass
			else:
				smut=smut+str(i+1)+seq[i]+" "               
                        
	return smut



######READING FASTA FILES and Translating + counting frequency#######
dProt={}
nFasta=0
for line in fFasta:
	if line.startswith('@'):
		nFasta+=1
	else:
		if len(ParL3)==0:
			monoseq = "".join(line.strip().split("\t"))
			protseq = translate(monoseq)
			if len(protseq)!=0:
				dProt = freq_dict(dProt,protseq)
		else:
			h3seq = translate(line.strip().split("\t")[0])
			l3seq = translate(line.strip().split("\t")[1])
			if (len(h3seq)!=0) and (len(l3seq)!=0):
				protseq = h3seq + "\t" + l3seq
				dProt = freq_dict(dProt,protseq)
fFasta.close()


#Write translated protein seq

fProtHtml.write("<!DOCTYPE html>\n<html>\n<head>\n<title>{0}</title>\n</head>\n<body>\n<p style='font-family:Courier New'>".format(dirpath))
fProtHtml.write("<strong>{}</strong><br><br>".format(ParH3+" "+ParL3))


nProt=0
nThrProt=0
nHeatProt=0
dmut={}
#in dProt dictionary, key is protein sequence and val is the frequency the read was observed. Eg: dProt["ASDED"]=24984
for key, val in sorted(dProt.items(), key=lambda item:item[1], reverse = True):
	nProt+=int(val)

	if dProt[key] >= nThr:
		nThrProt+=int(val)
		
		#smut has 1-based numbering of position + mutation + " " + so on. Eg. "51V 4A "
		#If L3 is present, the postion of mutation is position in L3 + len(H3) + 1
		smut = count_mutation(key,ParH3,ParL3)

		#Dictionary to write out Mutation in tab delimited format, only admitting one with single mutation as that was the library design
		#smut = position + amino acid + " "
		#dmut[amino acid] = [(position1(0-base numbering),frequency1),(position2,freqeuncy),...]
		if (len(smut.strip().split())==1) and (smut.find("_")==-1):
			nHeatProt += int(val)
			try:
				dmut[smut[-2:-1]].append((int(smut[:-2])-1,int(val)))
			except KeyError:
				dmut[smut[-2:-1]]=[(int(smut[:-2])-1,int(val))]				


		#This file (output in .prot format), output is protein sequence written number of times it was observed in the dataset
		#for j in range(int(dProt[key])):
		#	fProt.write("{0}\n".format(key))

		#Everything is tablimited in this output, .txt	
		tabmut=""
		for k in smut.strip().split(" "):
			tabmut+=k[:-1]+"\t"+k[-1:]+"\t"
		fProtTxt.write("{0}\t{1}\t{2}\t{3}\n".format(key,int(val),len(smut.strip().split(" ")),tabmut.strip()))
		
		#Writes html code for easy viewing with mutations colored red
		i=0
		for posmut in smut.split():
			pos,mut=int(posmut[:-1])-1,posmut[-1]	
			fProtHtml.write('{0}<span style="color:#FF0000">{1}</span>'.format(key[i:pos],key[pos]))
			i=pos+1
		fProtHtml.write('{0}&emsp;&emsp;{1}&emsp;&emsp;{2}<br>'.format(key[i:],dProt[key],smut))
#fProt.close()
fProtTxt.close()
fProtHtml.write("</p>\n</body>\n</html>")        
fProtHtml.close()

#Write  down whole parental seqeunce, tab separated
fProtHeat.write('WT')
for i in ParH3:
	fProtHeat.write('\t{}'.format(i))

if len(ParL3)!=0:
	fProtHeat.write('\t')
	for i in ParL3:
		fProtHeat.write('\t{}'.format(i))
fProtHeat.write('\n')


for i in sorted(dmut.keys()):
	fProtHeat.write(i)
	l=0
	for j,k in sorted(dmut[i]):
		normk=k*1000/nHeatProt
		if l ==0:
			blank="\t0"*(j+1)
			fProtHeat.write("{}{}".format(blank[:-1],normk))
		else:
			blank="\t0"*(j-l)
			fProtHeat.write("{}{}".format(blank[:-1],normk))
		l=j
	if len(ParL3)==0:
		blank = "\t0"*(len(ParH3)-l)
	else:
		blank = "\t0"*(len(ParH3+ParL3-1)-l)
	fProtHeat.write("{}\n".format(blank[:-1]))
		
	
fProtHeat.close()



print ("")
print ("{:60}{}".format("Number of Records in Fasta                             :",nFasta))
print ("{:60}{}".format("Number of Translated Records                           :",nProt))
print ("{:60}{}\t({})".format("Number of Records above Threshold and Threshold        :",nThrProt,nThr))
print ("{:60}{}".format("Number of Records with Single Mutation (HeatMap)       :",nHeatProt))



