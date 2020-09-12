import shapely
from shapely.geometry import shape, Point, mapping, MultiPoint
from shapely.ops import linemerge
import fiona
import networkx as nx
import csv
import scipy
from rtree import index
from math import exp, sqrt
from datetime import datetime
import itertools
from collections import Counter, defaultdict
import pandas as pd
import numpy as np
import itertools
from scipy.special import comb
import math
from numpy.linalg import matrix_power
import os
from os.path import isfile, join
import json
from pandas.io.json import json_normalize
from pyproj import Proj, transform
import utm
import geopandas as gpd
from shapely import wkt
import matplotlib.pyplot as plt
from shapely.geometry import box, Polygon
import osmnx as ox
import ast
import requests
from shapely.ops import nearest_points
import geopandas as gpd
from itertools import permutations
from lxml import objectify
import xml.etree.cElementTree as et
import KMPalgo as kmp

FMT = '%H:%M:%S'
def DTtoFloat(dstring):
    x = datetime.strptime(dstring, FMT)
    return x.hour * 3600 + x.minute * 60 + x.second

def LinkCounter(L,l): ## to check wheater l is a proper subset of L 
    indexMatch = list(kmp.KnuthMorrisPratt(L,l))
    if indexMatch:
        return 1#len(indexMatch)
    else:
        return 0

def ODCounter(L,lod): ## lod , a list of lists , contains multiple or single lists
    indexMatch = [] ## 
    for l in lod:## to check wheater l is a proper subset of L 
        indexMatch = indexMatch + list(kmp.KnuthMorrisPratt(L,l))
    
    if indexMatch:
        return len(indexMatch)
    else:
        return 0

## check wheather an OD contains other ODs
def ODcontainer(toSelectOne): #Pass a dataframe
    for item in range(len(toSelectOne)):
        print(item)
    return indexOf

AF = pd.read_csv('ODFinal.csv',index_col=0)
linksOD = [ast.literal_eval(item) for item in AF.index.to_list()]
AF = AF.loc[:, (AF != 0).any(axis=0)]
nameOD = [ast.literal_eval(item) for item in AF.dtypes.index.to_list()] #list(AF.dtypes.index[1::])

A = AF.values
NumColA = len(A[0])
NumRowA = len(A)


plt.rcParams.update({'font.size': 20})


errList = []
scr = pd.read_csv('ScatOsmPointRoutes.csv',index_col = 0)
scr.ScatRoute = scr.ScatRoute.apply(lambda x: ast.literal_eval(x))
scr.OsmRoute = scr.OsmRoute.apply(lambda x: ast.literal_eval(x))
scr['No'] =  scr.ScatRoute.apply(lambda x: len(x))
scr['NoLinks'] = scr.ScatRoute.apply(lambda x: len(x[0])) #scr.apply(lambda x: len(x.ScatRoute[0]),axis=1)
scr['MaxNoLinks'] = scr.ScatRoute.apply(lambda x: max([len(i) for i in x])) ## maximum links involved
scr['LongestRoute'] = scr.ScatRoute.apply(lambda x: x[np.argmax([len(i) for i in x])])

scr.index = list(zip(scr.O.to_list(),scr.D.to_list()))



path = "Smarts Random OD/"
abspath = os.path.abspath(path)
filenames = os.listdir(path)
dfcols = ['Vid','type','start','osmid']
simdf = pd.DataFrame(columns=dfcols)
InitTime = 0
for file in filenames:
    cfile = path + file
    xml = et.parse(cfile)
    for i in xml.iter(tag='vehicle'):
        osmids = []
        vid = i.items()
        tempListNodes = list(i.iter(tag='node'))
        for o in tempListNodes:
            osmids.append(o.items()[0][1])
        simdf = simdf.append(pd.Series([vid[0][1],vid[1][1],InitTime + float(vid[2][1]),osmids], index=dfcols),ignore_index=True)
    InitTime = InitTime + 950
simdf = simdf[simdf.type == 'CAR']

simdf["O"] = simdf["Vid"].apply(lambda x: int(x[0:4]))
simdf["D"] = simdf["Vid"].apply(lambda x: int(x[5:9]))
#simdf["hod"] = simdf["Vid"].apply(lambda x: x[10:12])
#simdf["Q"] = simdf["Vid"].apply(lambda x: x[13:15])
simdf["OD"] = list(zip(simdf.O.tolist(),simdf.D.tolist()))
#simdf["time_interval"] = simdf["Vid"].apply(lambda x: x[10:18])
#create a linkdf from scr that contains

linkdf = scr[scr.NoLinks == 2].copy()
linkdf.drop(['No','NoLinks'],axis=1,inplace=True)

## getting link counts from df

minTime = 0
tempTime = 0
maxTime = simdf.iloc[-1].start
intervals = 300 ## create interval of 300 seconds
## for each interval, create link count

while minTime < maxTime:
    tempTime = minTime + 900
    simdfTrim = simdf[(minTime <= simdf.start) & (simdf.start <= tempTime)]
    minTime = tempTime
    allL = simdfTrim.osmid.to_list()
    linkdf[str(tempTime)] = [0] * len(linkdf)
    scr[str(tempTime)] = [0] * len(scr)

    print(len(allL))
    print('-----S T A R T ------')
    for idex,temp in enumerate(allL):
        #print(idex)
        L = list(map(int, temp))
        linkdf[str(tempTime)] = linkdf[str(tempTime)] + linkdf.OsmRoute.apply(lambda x: LinkCounter(L,x[0]))
        redundent = scr.OsmRoute.apply(lambda x: ODCounter(L,x)) ## by doing this you r redundantly adding ODs for sub Ods
        ## hence you have consider the Only OD with laregest possible route
        idxs = redundent[(redundent > 0)].index.to_list()
        if idxs:
            toSelectOne = scr.loc[idxs]
            OnlyIndexToUpdate = toSelectOne['MaxNoLinks'].idxmax()
            value = redundent.get(OnlyIndexToUpdate)
            ## select the OD that covers other redundent ODs
            scr.at[OnlyIndexToUpdate,str(tempTime)] = scr.at[OnlyIndexToUpdate,str(tempTime)] + 1 

columnsWithLinkcounts = linkdf.columns.to_list()[6::]
Y_all_time = linkdf[columnsWithLinkcounts].values.tolist()
#-------------------- Calculate inter arrival rate from linkdf
if 1 == 1:
    Lam = (np.ones(NumColA))
    LamAll = np.ones(NumColA)

    Y = np.asarray(Y_all_time)  ## put link counts here Y= K * Number of link matrix
    CovY = np.cov(Y,bias=True)

    #np.where(~A.any(axis=0))[0]
    Y_bar = np.vstack((np.mean(Y,axis=1)))
    S_all = np.hstack((np.triu(CovY)))
    S = np.vstack((S_all[S_all!=0])) 

    tempB0 = np.multiply(A[0::], A[0])
    for CountB in range (1,len(A)):
        tempB0 = np.vstack((tempB0,np.multiply(A[CountB::], A[CountB])))
        #print(len(tempB0))
    B = tempB0 * 0.2 # regularization
    S = S * 0.2 # regularization
    n_B = np.delete(B, np.where(~B.any(axis=1)),0)
    n_S = np.delete(S, np.where(~B.any(axis=1)),0)

    Pos_S = n_S
    Pos_S[Pos_S < 0] = 0

    F_B = np.delete(n_B, np.where(~Pos_S.any(axis=1)),0)
    F_S = np.delete(n_S, np.where(~Pos_S.any(axis=1)),0)

    #AFinal = np.vstack((A,F_B))
    #YFinal = np.vstack((Y_bar,F_S))
   
    AFinal = A
    YFinal = Y_bar
    TotalRow = len(YFinal)
    TempLam = np.ones(NumColA)
    #Y[0:4:3, 0:4:3]
    for numIter in range (0,50):
        print(numIter)
        for Countj in range (0,NumColA):
            sTerm = np.zeros(TotalRow)
            for row in range (0,TotalRow):
                sTermD = np.sum(np.multiply(AFinal[row,:],Lam))
                #print(sTermD)
                if sTermD == 0.0:
                    print(row)
                    break
                sTermN = float(AFinal[row,Countj]) * float(YFinal[row])
                #print(sTermN/sTermD)
                sTerm[row] = sTermN/sTermD
            #print(np.sum(sTerm))
            TempLam[Countj] = float(Lam[Countj]) * float((np.sum(sTerm))/np.sum(AFinal[:,Countj]))
            #print("----E N D O F L O O P --------")
            #print(Countj)
        Lam = TempLam
        #LamAll = np.vstack((LamAll,Lam))
    iar = TempLam

## ------------------ compare with scr ------------ ##
if 1 == 1:
    columnsWithODflows = scr.columns.to_list()[8:-1]
    Xod = scr[columnsWithODflows].copy()
    Xod['MeanOD'] = Xod.mean(axis=1)
    Xod['Estimated'] = iar
    Xod['Estimated'] = Xod.Estimated.apply(lambda x: round(x,2) if (x > 0.0) else 0)
    Xod['O'] = scr.O
    Xod['D'] = scr.D
    xod = Xod[Xod.MeanOD >= 1]
    #xod = Xod[(Xod.O != 4389) & (Xod.O != 2905)]
    #xod = xod[(xod.O != 4403) & (xod.D != 4403)]
    #xod = xod[(xod.O != 2904) & (xod.O != 4403)]
    xob = xod.MeanOD.values
    xeq = xod.Estimated.values
    rmse = ((xod.MeanOD - xod.Estimated) ** 2).mean() ** .5
    mape = abs(np.mean((xod.MeanOD - xod.Estimated) / xod.MeanOD * 100))
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    #x = np.arange(len(sample))
    ax.scatter(xeq,xob)
    x = np.linspace(0, max(xob), 5) # constructs a numpy array of [0.0, 1.0, ... 10.0]
    ax.plot(x, x, linestyle='-',c='r', label='Slope 1 line')
    ax.legend(loc='lower right')
    ax.set_ylabel('OD flows observed')
    ax.set_xlabel('OD flows estimated')
    #Title = 'Estimated travel time comparison with Sygic data (Q1 of year 2016)'
    #ax.set_title(Title)
    #ax.set_ylim(-1,15)
    plt.show()
    Xod.to_csv('Xod.csv',header=True)

if 1 == 1:
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    x = np.arange(0,len(xod))
    xod = xod.sort_values(by=['MeanOD'])
    xob = xod.MeanOD.to_list()
    xeq = xod.Estimated.to_list()
    ax.plot(x,xeq,ls='-',marker='o',c='r',label='Estimated OD flow')
    ax.plot(x,xob,ls='-',marker='*',c='b',label='Observed OD flow')
    ax.legend(loc='upper left')
    ax.set_ylabel('Number of vehicles travelled')
    ax.set_xlabel('OD pair')
    #ax.set_xticks(xod.index.to_list())
    plt.show()
