#!/usr/bin/python3

import os
import sys
import numpy as np
import pandas as pd


#Variables here
#################################
##YOU NEED TO EDIT DIRPATH#######
#################################
dirpath = sys.argv[1]

#update Parental seq!
sParSeq = sys.argv[2].replace('_',' ')
lParSeq = sParSeq.split()

#Takes in mutation positions allowed by NNK design
sMutPos = sys.argv[3]

#List of Amino acids in alphabetical order
lAmino = sorted(['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V'])

#Creating Array filled with value of 0 to save mutation frequency information in array.
arrMutCount=np.zeros([len(lAmino),len(sParSeq)])

#Reads with frequency <= nThr will be ignored in the output
nThr = int(sys.argv[4])

#Number of NNK mutation
nNNK = int(sys.argv[5])

#opens files to read and write
for fName in os.listdir(dirpath):
	if fName.endswith(".seq"):
		outName=fName[:-4]
		fSEQ = open(os.path.join(dirpath,fName),'r')


#Html file to write
fHtml = open(os.path.join(dirpath,outName+"_"+".html"),'w')

fHtml.write("<!DOCTYPE html>\n<html>\n<head>\n<title>{0}</title>\n</head>\n<body>\n<p style='font-family:Courier New'>".format(dirpath))
fHtml.write("<strong>{}</strong><br><br>".format(' '.join(lParSeq)))



#####################
##Reading .SEQ file and saving output in dictionary of arrays
nInit,nUnique,nUniqueFin,nFin,nInsDel,nWeird,nlowerThr,nWT=0,0,0,0,0,0,0,0
dSeq={}
for line in fSEQ:
	lLine=line.split('\t')
	sLine=' '.join(lLine[:-1])	
	nCount=int(lLine[-1])
	f=0	#Flage indicator to break when record is not needed

	nInit += nCount
	nUnique +=1
	
	#Cheks if the count is above threshold
	if nCount < nThr:
		nlowerThr+=nCount
		f=1
		continue
	if f==1: continue	
	#Check if Sequence have same length as ParSeq. For NNK lib, theoretically no insert/deletion should be observed
	for i in range(len(lParSeq)):
		par,seq=lParSeq[i], lLine[i]
		if len(par) != len(seq):
			nInsDel += nCount
			f=1
			break
	if f==1: continue
	#Checks if there is right number, position (0-base) of mutation allowed by NNK primer design
	lMut=[m for m, (s1,s2) in enumerate(zip(sParSeq,sLine)) if s1 != s2]
	if (len(lMut) > nNNK) and (len(lMut) !=0):
		nWeird+= nCount
		f=1
		continue
	for m in lMut:
		if sMutPos[m] != 'M':
			nWeird += nCount
			f=1
			break
	if f==1: continue

	#Count for PARENTAL. If you do not want to count these, indent these out.
	if sParSeq==sLine:
		f=1
		nWT += nCount
		nFin += nCount
		nUniqueFin +=1
		dSeq[sLine]=nCount
		for i,j in enumerate(sParSeq):
			if j == " ": 
				pass
			else:
				nAmino=lAmino.index(j)
				arrMutCount[nAmino,i] += nCount
		fHtml.write('{}&emsp;&emsp;{}<br>'.format(sLine,nCount))
	
	#Check if flagged
	if f==1: continue

	#Save in arrMutCount Array.
	#Save smut to write to html
	sMut=""
	nFin += nCount
	nUniqueFin +=1
	i=0
	dSeq[sLine]=nCount
	for m in lMut:
		nAmino = lAmino.index(sLine[m])
		arrMutCount[nAmino,m] += nCount
		fHtml.write('{}<span style="color:#FF0000">{}</span>'.format(sLine[i:m],sLine[m]))
		i=m+1
		sMut += str(m+1)+sLine[m]+" "
	fHtml.write('{}&emsp;{}&emsp;&emsp;{}'.format(sLine[i:],nCount,sMut[:-1]))
	fHtml.write("<br>")
	for m2,m3 in enumerate(sLine):
		if m2 not in lMut:
			nAmino=lAmino.index(sParSeq[m2])
			arrMutCount[nAmino,m2] += nCount

fSEQ.close()
fHtml.close()


#Saving Arrays in CVS. adding 1 for mathematical reason.
#MutFreq is frequency of each amino acid at each position
arrMutFreq=(arrMutCount+1)*100/nFin

dfMutCount = pd.DataFrame(arrMutCount, index=lAmino, columns=list(sParSeq))
dfMutFreq = pd.DataFrame(arrMutFreq, index=lAmino, columns=list(sParSeq))

dfMutCount.to_csv('%s/%s_count.csv'%(dirpath, outName))
dfMutFreq.to_csv('%s_freq.csv'%(outName))

#Seq and freqeuncy. This output looks at sequence as whole
lSeqFreq=[]
for i in dSeq.values():
	f=i*100/nFin
	lSeqFreq.append(f)
dfSeqFreq = pd.DataFrame(lSeqFreq, index=dSeq.keys(), columns=["Freq"])
dfSeqFreq.to_csv('%s_seqfreq.csv'%outName)

			
print (outName)
print ("{:70}{}{}".format("Number of Records in SEQ",": ",nInit))
print ("{:70}{}{}".format("Number of Ins/Del Mutation",": ",nInsDel))
print ("{:70}{}{}\t({})".format("Number of Records Lower Than Threshold (Threshold)",": ",nlowerThr,nThr))
print ("{:70}{}{}".format("Number of Records with diff number/position mut from designated",": ",nWeird))
print ("{:70}{}{}".format("Number of WT Records",": ",nWT))
print ("{:70}{}{}".format("Number of Final Records",": ",nFin))
print ("{:70}{}{}".format("Number of Initial Unique Records",": ",nUnique))
print ("{:70}{}{}\n".format("Number of Final Unique Records",": ",nUniqueFin))


