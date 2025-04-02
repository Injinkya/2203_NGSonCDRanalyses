#!/usr/bin/python3

import os
import sys
import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

dDF={}
for fName in os.listdir():
	if fName.endswith("seqfreq.csv"):
		dDF[fName[:4]] = pd.read_csv(r'%s'%(fName), index_col=0, header =0)

def get_common(lib1,lib2):
	#dfmerge = pd.merge(dDF[lib1],dDF[lib2], left_index=True, right_index=True)
	dfconcat = pd.concat([dDF[lib1],dDF[lib2]],axis=1)
	#dfmerge.columns =[lib1,lib2]
	dfconcat.columns =[lib1,lib2]
	#dfmerge.to_csv('%s_%s_merge.csv'%(lib1,lib2))
	dfconcat.to_csv('%s_%s_concat.csv'%(lib1,lib2))
	return  dfconcat #, dfmerge



dfcon12 = get_common("GK-1","GK-2")
dfcon35 = get_common("GK-3","GK-5")
dfcon46 = get_common("GK-4","GK-6")

dfcon1235 = pd.concat([dfcon12,dfcon35],axis=1)
dfcon1246 = pd.concat([dfcon12,dfcon46],axis=1)
dfcon123546 = pd.concat([dfcon1235,dfcon46],axis=1)

dfcon1235.to_csv("GK-1235_concat.csv")
dfcon1246.to_csv("GK-1246_concat.csv")
dfcon123546.to_csv("GK-123546_concat.csv")


df46not12 = dfcon46[~dfcon46.index.isin(dfcon12.index)]
df35not12 = dfcon35[~dfcon35.index.isin(dfcon12.index)]
df46not12.to_csv("GK46binder.csv")
df35not12.to_csv("GK35binder.csv")
df46unique = df46not12[~df46not12.index.isin(df35not12.index)]
df35unique = df35not12[~df35not12.index.isin(df46not12.index)]
df46unique.to_csv("GK46binder_unique.csv")
df35unique.to_csv("GK35binder_unique.csv")

	
