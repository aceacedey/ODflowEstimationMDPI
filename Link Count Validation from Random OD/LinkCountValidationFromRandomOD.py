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
    ## check for consecutive duplicates in l
    #nwl = [i[0] for i in itertools.groupby(l)]
    #print(nwl)
    indexMatch = list(kmp.KnuthMorrisPratt(L,l))
    if indexMatch:
        return 1#len(indexMatch)
    else:
        return 0

def ODCounter(L,lod): ## lod , a list of lists , contains multiple or single lists
    indexMatch = [] ## 
    for l in lod:## to check wheater l is a proper subset of L 
    ## check for consecutive duplicates in l
        #nwl = [i[0] for i in itertools.groupby(l)]
        #print(nwl)
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
scr['SimOD'] = [0] * len(scr)

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


allL = simdf.osmid.to_list()
#create a linkdf from scr that contains
if 1 == 1:
    linkdf = scr[scr.NoLinks == 2].copy()
    linkdf.drop(['No','NoLinks'],axis=1,inplace=True)
    linkdf['SimLC'] = [0] * len(linkdf)
    ## getting link counts from df
    print(len(allL))
    print('-----S T A R T ------')
    for idex,temp in enumerate(allL):
        #print(idex)
        L = list(map(int, temp))
        linkdf['SimLC'] = linkdf['SimLC'] + linkdf.OsmRoute.apply(lambda x: LinkCounter(L,x[0]))
        redundent = scr.OsmRoute.apply(lambda x: ODCounter(L,x)) ## by doing this you r redundantly adding ODs for sub Ods

        ## hence you have consider the Only OD with laregest possible route
        idxs = redundent[(redundent > 0)].index.to_list()
        
        if idxs:
            toSelectOne = scr.loc[idxs]
            OnlyIndexToUpdate = toSelectOne['MaxNoLinks'].idxmax()
            value = redundent.get(OnlyIndexToUpdate)
            ## select the OD that covers other redundent ODs
            #scr['SimOD'][scr.index == OnlyIndexToUpdate] = scr['SimOD'][scr.index == OnlyIndexToUpdate] + 1
            scr.at[OnlyIndexToUpdate,'SimOD'] = scr.at[OnlyIndexToUpdate,'SimOD'] + 1 
            
    print('-----E N D ------')
    Y_ob = linkdf.SimLC.to_list()
    X = scr.SimOD.to_list()
    Y_eq = np.matmul(A,X)

  
    resdf = pd.DataFrame({'equ':list(Y_eq),'ob':list(Y_ob),'dif':list(Y_eq - Y_ob)},index=linkdf.index)
    fdf = resdf.copy()
linkdf['Estimated'] = list(Y_eq)
linkdf.to_csv('Ylc.csv')
##if 1 == 1:
##    b = [(4389,4403),(4389,4523),(2904,2935),(2919,4405),(2904, 2905),(2904, 2913),(4404, 4403)]
##    fdf = fdf[~fdf.index.isin(b)]
##    yeq = fdf.equ.values
##    yob = fdf.ob.values
##    r = np.sqrt(np.mean(yeq-yob)**2)
##    #plt.scatter(yeq,yob)
##    
##    #plt.show()
##if 1 == 1: # plot link validation: 
##
##    fig = plt.figure()
##    ax = fig.add_subplot(1, 1, 1)
##    #x = np.arange(len(sample))
##    ax.scatter(yeq,yob)
##    x = np.linspace(0, max(yob), 5) # constructs a numpy array of [0.0, 1.0, ... 10.0]
##    ax.plot(x, x, linestyle='-',c='r', label='Slope 1 line')
##    ax.legend(loc='lower right')
##    ax.set_ylabel('Link counts observed')
##    ax.set_xlabel('Link counts from estimated OD flow')
##    #Title = 'Estimated travel time comparison with Sygic data (Q1 of year 2016)'
##    #ax.set_title(Title)
##    #ax.set_ylim(-1,15)
##    plt.show()
##    fdf.to_csv('1000fold.csv',header=True)

if 1 == 1:
   # resdf = pd.read_csv("
   # yod = linkdf[(linkdf.O != 4389) & (linkdf.O != 2905)]
    #b = [(4389,4403),(4389,4523),(2904,2935),(2919,4405),(2904, 2905),(2904, 2913),(4404, 4403)]
    #yod = linkdf[~linkdf.index.isin(b)][['SimLC','Estimated']]#/4
    yod = linkdf[['SimLC','Estimated']]
    yod = yod[(yod.SimLC != 0) & (yod.Estimated != 0)]
    #xod = xod[(xod.O != 4403) & (xod.D != 4403)]
    #yod = yod[(yod.O != 2904) & (yod.O != 4403) ]
    yob = yod.SimLC.values
    
    yeq = yod.Estimated.values
    rmse = ((yod.SimLC - yod.Estimated) ** 2).mean() ** .5
    mape = abs(np.mean((yod.SimLC - yod.Estimated) / yod.SimLC)) * 100
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    #x = np.arange(len(sample))
    ax.scatter(yeq,yob,marker='*')
    x = np.linspace(0, max(yob), 5) # constructs a numpy array of [0.0, 1.0, ... 10.0]
    ax.plot(x, x, linestyle='-',c='r', label='Slope 1 line')
    ax.legend(loc='lower right')
    ax.set_ylabel('Link counts observed')
    ax.set_xlabel('Link counts estimated')
    #Title = 'Estimated travel time comparison with Sygic data (Q1 of year 2016)'
    #ax.set_title(Title)
    #ax.set_ylim(-1,15)
    plt.show()
    #Xod.to_csv('Xod.csv',header=True)
