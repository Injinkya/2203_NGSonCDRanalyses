#!/usr/bin/python3

import os
import sys
from Bio import SeqIO 

#Variables here
#################################
##YOU NEED TO EDIT DIRPATH#######
#################################
dirpath = sys.argv[1]


#Update the flanking SEQ
FW = sys.argv[2]
RV = sys.argv[3]
#FW="mono"
#RV="mono"

#H1 SGFT_H1_RQAP	GLEW_H2_RFTI
if FW =="H1H2":
	NterFW="GGCTTCACA"
	CterFW="CGTCAGGCC"
	NterFW1="TGGAATGG"
	CterFW1="CGTTTCACT"

#H3 YCAR_H3_WGQG
elif FW =="H3":
	NterFW="TGTGCTCGC"
	CterFW="TGGGGTCAA"

#L3 FATYYC_L3_FGQGT
elif FW =="L3":
	NterFW="TATTACTGT"
	CterFW="TTCGGACAG"
elif FW =="mono":
	RV="mono"
	NterFW="CCGACTAGC"
	CterFW="ATTAGCGGC"
else:
	print ("FW has to be either 'H1H2' or 'H3' or 'L3'")
	sys.exit()

if RV =="H3":
	NterRV= "TGTGCTCGC"
	CterRV= "TGGGGTCAA"
elif RV =="L3":
	NterRV= "TATTACTGT"
	CterRV= "TTCGGACAG"
elif RV =="mono":
	NterRV= "CCGCGACC"
	CterRV= "TACCGTAC"
else:
	print ("RV has to be either 'H3' or 'L3' or 'mono'")
	sys.exit()

#opens files to read and write
for fName in os.listdir(dirpath):
	if fName.endswith('.fastq'):
		if "_R1_" in fName:
			forFQ = open(os.path.join(dirpath,fName),'r')		
			a=fName.find("_R1_")
			outName=(fName[:a])+".seq"
			outSEQ = open(os.path.join(dirpath,outName),'w')
		elif "_R2_" in fName:
			revFQ = open(os.path.join(dirpath,fName),'r')


#function to reverse complement the DNA seq
def rev_compl(seq):
	compl={'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
	return ''.join(compl[i] for i in seq[::-1])

#protein translate
def translate(seq):

	table= {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M','ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K','AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L','CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
                'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q','CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V','GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
                'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E','GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
                'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S','TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
                'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_','TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',}   
	protein =""
	if len(seq)%3 == 0:
		for i in range(0, len(seq), 3): 
			codon = seq[i:i + 3]
			if codon.find("N")==-1:
				protein+= table[codon]
			else:
				return ""
	else:
		pass

	return protein




######READING FASTQ FILES#########
#reading forward file in to dfor dictionary, trimming the seq to within the given NterFW, CterFW seq
dfor={}
dforno={}
nforFQ=0
nforFQno=0
nforFQqc=0
for record in SeqIO.parse(forFQ,'fastq'):
	nforFQ+=1
	name=record.id
	s=record.seq
	ps=record.letter_annotations["phred_quality"]
	try:
		d=s.index(NterFW)+len(NterFW)
		e=s.index(CterFW)
		#e=d+150
		if min(ps[d:e]) >=1:
			dfor[name]=s[d:e]
		else:
			nforFQqc+=1
			continue
	except ValueError:
		nforFQno+=1
		dforno[name]=s
forFQ.close()

#reading reverse file, do reverse complment and save in drev dictionary, trimming the seq to NterRV and CterRV
drev={}
drevno={}
nrevFQ=0
nrevFQno=0
nrevFQqc=0
for record in SeqIO.parse(revFQ,'fastq'):
	nrevFQ+=1
	name=record.id
	s=rev_compl(record.seq)
	ps=record.letter_annotations["phred_quality"]
	try:
		d=s.index(NterRV)+len(NterRV)
		e=s.index(CterRV)
		#d=e-150
		if min(ps[d:e]) >= 1:
			drev[name]=s[d:e]
		else:
			nrevFQqc+=1
			continue
	except ValueError:
		nrevFQno+=1
		drevno[name]=s
revFQ.close()

#write files
#Unique ID for each sample \n Forward read \t Reverse read \n
nPaired,nNotPaired,nTranslated,nStop=0,0,0,0
dCount={}
for key in dfor:
	try:
		seq1,seq2=str(dfor[key]),str(drev[key])

		#nseq2=seq2.index(seq1[-7:])+7
		nPaired +=1
		
		protseq1=translate(seq1)
		protseq2=translate(seq2)
		#Check that both forward, reverse seq are translated and save only those.
		if len(protseq1)!=0 and len(protseq2)!=0:
			nTranslated +=1
			protseq=protseq1+protseq2
		else:
			continue
		
		if protseq.find("_")==-1:
			try:	
				dCount[protseq]+=1			
			except KeyError:
				dCount[protseq]=1
		else:
			nStop+=1

	except KeyError:
		nNotPaired+=1
	except ValueError:
		nNotPaired+=1


###Writing output file
for protseq,freq in sorted(dCount.items(), key=lambda item:item[1], reverse=True):
	outSEQ.write("{}\t{}\n".format(protseq,freq))
outSEQ.close()





##these files are just to check rejected files without flanking sequences.
# Delete # if you want these files written
#outFQforno = open(os.path.join(dirpath,fName[:a]+"for_no.fa"),'w')
#outFQrevno = open(os.path.join(dirpath,fName[:a]+"rev_no.fa"),'w')
#for key in dforno:
#	outFQforno.write("{0}\n{1}\n".format(key,dforno[key]))
#for key in drevno:
#	outFQrevno.write("{0}\n{1}\n".format(key,drevno[key]))
#outFQforno.close()
#outFQrevno.close()



print ("{}".format(dirpath[:20]))
print ("{:60}{}{}".format("Number of Fastq Forward Records",": ",nforFQ))
print ("{:60}{}{}".format("Number of Fastq Reverse Records",": ",nrevFQ))
print ("{:60}{}{}".format("Number of Forwards wo adaptor seq",": ",nforFQno))
print ("{:60}{}{}".format("Number of Reverses wo adaptor seq",": ",nrevFQno))
print ("{:60}{}{}".format("Number of Forwards one/more base with <20 Phred score",": ",nforFQqc))
print ("{:60}{}{}".format("Number of Reverses one/more base with <20 Phred score",": ",nrevFQqc))
print ("{:60}{}{}".format("Number of Not Paired Records",": ",nNotPaired))
print ("{:60}{}{}".format("Number of Paired Records",": ",nPaired))
print ("{:60}{}{}".format("Number of Stop Mutations",": ",nStop))
print ("{:60}{}{}".format("Number of Translated records",": ",nTranslated))
print ("{:60}{}{}\n".format("Number of Unique Records",": ",len(dCount.keys())))


