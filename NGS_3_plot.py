#!/usr/bin/python3

import os
import sys
import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

#Variables here
#################################
##YOU NEED TO EDIT DIRPATH#######
#################################
sOri=sys.argv[1]
sLib=sys.argv[2]

#update Parental seq!
sParSeq = sys.argv[3].replace('_',' ')
lParSeq = sParSeq.split()

typePlot= sys.argv[4]

#Normalizing Cell!
norm=(1,-1)


lCDRadj=[[0,0],[0,0],[0,0]]
lCDRpos =[]
j=0
for i, seq in enumerate(lParSeq):
	start=j+lCDRadj[i][0]
	end=j+len(seq)+lCDRadj[i][1]
	lCDRpos.append((start,end))
	j=j+len(seq)

enmin, enmax = -2, 2

#List of Amino acids in alphabetical order
lAmino = sorted(['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V'])

#opens files to read and write
dArr={}
dEMM={}

#Reading all freq_threshold.csv file into numpy and storing in dictionary format.
# {....freq_1 : numpy}
for fName in os.listdir():
	if fName.endswith("freq.csv"):
		outName=fName[:-4]		
		arr=pd.read_csv(r'%s'%(fName),index_col=0,header=0).to_numpy()
		j=0	
		dArr[outName]=np.empty((20,0))
		for i,seq in enumerate(lParSeq):
			l=len(seq)+j
			dArr[outName]=np.append(dArr[outName], arr[:,j:l], axis=1)
			j=l+1
		if sOri in outName:
			arrOri=dArr[outName]

#Break down each array to each CDR and save both freq and enrichment
dArrFin={}
l=[]
#ll=[]
for name in dArr.keys():
	#Enrichment is calculated here
	if typePlot =="log10":
		arren=np.log10(dArr[name]/arrOri)
		arrfr=np.log10(dArr[name])
	elif typePlot =="log2":
		arren=np.log2(dArr[name]/arrOri)
		arrfr=np.log2(dArr[name])
	normen = arren[norm[0],norm[1]] 
	normfr = arrfr[norm[0],norm[1]]
	
	#adjusting data for enrichment within bonundaries
	arren=np.where(arren>=enmax, enmax, arren-normen)
	arren=np.where(arren<=enmin, enmin, arren)
	

	#no need to adjust for frequency, so letting it be
	arrfr=arrfr-normfr

	for j, (start,end) in enumerate(lCDRpos):
		fname=name+"_%s"%j
		ename=fname.replace("freq","enrich")
		dArrFin[fname]=arrfr[:,start:end]
		dArrFin[ename]=arren[:,start:end]
	if sLib in name:
		l.extend([np.max(arrfr),np.min(arrfr)])
	#	ll.extend([np.max(arren),np.min(arren)])
#evmin,evmax,fvmin,fvmax=min(ll),max(ll),min(l),max(l)
fvmin,fvmax=min(l),max(l)


##Storing width for figure
ll=[]
for i in range(len(lParSeq)):
	ll.append(len(lParSeq[i]))
	


#Set the color for the heatmap
cmap_en = sb.diverging_palette(220, 20, as_cmap=True)
cmap_fr = "Oranges"

def plotHM(cmap_type, freq_enrich,mini, maxi):
	fig, ax = plt.subplots(1,ncols=len(lParSeq), figsize=(20,6), gridspec_kw={'width_ratios':ll, 'wspace':0.2})
	cbar_ax = fig.add_axes([.87,.2,.01,.55])
	cbar_ax.set_title(typePlot,pad=10)
	fig.subplots_adjust(bottom=0.1, top=0.9, left=0.07, right=0.8)
	for i in dArrFin.keys():
		if sLib in i and freq_enrich in i:
			j=int(i[-1])
			try:
				axx=ax[j]
			except TypeError:
				axx=ax

			arr = dArrFin[i]
			
			seq="".join(lParSeq)		

			plt.figure()	#(figsize = (len(lParSeq[j])/4,len(lAmino)/4))
			collist=list(seq)[lCDRpos[j][0]:lCDRpos[j][1]]
			df = pd.DataFrame(arr, index=lAmino, columns=collist)
			sm = sb.heatmap(df, linewidths=.4, cmap=cmap_type, square=True, ax=axx, vmin=mini, vmax=maxi, cbar_ax=cbar_ax)
			axx.set_xticklabels(collist, ha='center')
			axx.xaxis.tick_top() 
			axx.set_yticklabels(lAmino, rotation=0, ha='center')
			axx.tick_params(axis='y', pad=10)
			cbar_ax.set_title(typePlot, pad=10)
			sm.tick_params(length=0)
	fig.savefig('%s_%s.png'%(sLib,freq_enrich), dpi=300)
	


plotHM(cmap_fr, "freq", fvmin, fvmax)
plotHM(cmap_en, "enrich",enmin, enmax)


#for i in dArrFin.keys():
#	if sLib in i and "enrich" in i:
#		j=int(i[-1])
#		try:
#			axx=ax[j]
#		except TypeError:
#			axx=ax
#
#		arr = dArrFin[i]
#		plt.figure()	#(figsize = (len(lParSeq[j])/4,len(lAmino)/4))
#		collist=list(sParSeq.replace(" ","")[lCDRpos[j][0]:lCDRpos[j][1]])	
#		df = pd.DataFrame(arr, index=lAmino, columns=collist)
#		sm = sb.heatmap(df, linewidths=.4, cmap=cmap_en,square=True, ax=axx, vmin=-1, vmax=1, cbar_ax=cbar_ax)
#		axx.set_xticklabels(collist, ha='center')
#		axx.xaxis.tick_top() 
#		axx.set_yticklabels(lAmino, rotation=0, ha='center')
#		axx.tick_params(axis='y', pad=10)
#			
#		sm.tick_params(length=0)
#fig.savefig('%s_enrich.png'%sLib, dpi=300)
	




