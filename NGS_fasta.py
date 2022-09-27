#!/usr/bin/python3

import os
import sys


#Variables here
#################################
##YOU NEED TO EDIT DIRPATH#######
#################################
dirpath = sys.argv[1]


#Update the flanking SEQ
NterFW = sys.argv[2]
CterFW = sys.argv[3]
NterRV = sys.argv[4]
CterRV = sys.argv[5]


#opens files to read and write
for fName in os.listdir(dirpath):
	if fName.endswith('.fastq'):
		if "_R1_" in fName:
			forFQ = open(os.path.join(dirpath,fName),'r')		
			a=fName.find("_R1_")
			outName=(fName[:a])+".fa"
			outFQ = open(os.path.join(dirpath,outName),'w')
			outFQforno = open(os.path.join(dirpath,fName[:a]+"for_no.fa"),'w')
			outFQrevno = open(os.path.join(dirpath,fName[:a]+"rev_no.fa"),'w')
		elif "_R2_" in fName:
			revFQ = open(os.path.join(dirpath,fName),'r')


#function to reverse complement the DNA seq
def rev_compl(seq):
	compl={'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
	return ''.join(compl[i] for i in seq[::-1])


######READING FASTQ FILES#########
#reading forward file in to dfor dictionary, trimming the seq to within the given NterFW, CterFW seq
dfor={}
dforno={}
j=0
nforFQ=0
nforFQno=0
for line in forFQ:
	if line.startswith('@') and (j==0 or j >= 4):
		j=0
		b=line.strip().split(" ")[0]
		nforFQ+=1
	elif j == 1:
		c=line.strip()
		try:
			d=c.index(NterFW)+len(NterFW)
			e=c.index(CterFW)
			dfor[b]=c[d:e]
		except ValueError:
			nforFQno+=1
			dforno[b]=c
	j+=1
forFQ.close()


#reading reverse file, do reverse complment and save in drev dictionary, trimming the seq to NterRV and CterRV
drev={}
drevno={}
j=0
nrevFQ=0
nrevFQno=0
for line in revFQ:
	if line.startswith('@') and (j==0 or j >=4):
		j=0
		b=line.strip().split(" ")[0]
		nrevFQ+=1
	elif j==1:
		c=rev_compl(line.strip())
		try:
			d=c.index(NterRV)+len(NterRV)
			e=c.index(CterRV)
			drev[b]=c[d:e]
		except ValueError:
			nrevFQno+=1
			drevno[b]=c
	j+=1
revFQ.close()


#write files
#Unique ID for each sample \n Forward read \t Reverse read \n
nGoodRec,nBadRec=0,0
for key in dfor:
	try:
		outFQ.write("{0}\n{1}\t{2}\n".format(key,dfor[key],drev[key]))
		nGoodRec +=1
	except KeyError:
		nBadRec+=1
outFQ.close()

for key in dforno:
	outFQforno.write("{0}\n{1}\n".format(key,dforno[key]))
for key in drevno:
	outFQrevno.write("{0}\n{1}\n".format(key,drevno[key]))
outFQforno.close()
outFQrevno.close()




print ("{:50}{}".format("Number of Fastq Forward Records        :",nforFQ))
print ("{:50}{}".format("Number of Fastq Reverse Records        :",nrevFQ))
print ("{:50}{}".format("Number of Forwards wo givne seq        :",nforFQno))
print ("{:50}{}".format("Number of Reverses wo given seq        :",nrevFQno))
print ("{:50}{}".format("Number of Good Records                 :",nGoodRec))
print ("{:50}{}".format("Number of Records not paired           :",nBadRec))
