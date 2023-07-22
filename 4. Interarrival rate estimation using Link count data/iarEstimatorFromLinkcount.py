import shapely
from shapely.geometry import shape, Point, mapping, MultiPoint
from shapely.ops import linemerge
import fiona
import networkx as nx
import csv
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
from os import listdir
from os.path import isfile, join
import json
from pandas.io.json import json_normalize
import datetime
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
plt.rcParams.update({'font.size': 18})

def PlotTraj(traj,pointsp,ogp):
    lon = []
    lat = []
    for item in pointsp:
        lon.append(item[0])
        lat.append(item[1])
    fig, ax = ox.plot_graph(ogp, show=False, close=False)
    
    if traj:
        x1 = np.arange(0,len(traj))
    
        osmN,osmE = ox.graph_to_gdfs(ogp,edges=True,nodes=True)
        poi = [(osmN.x[i],osmN.y[i]) for i in traj]
        poix = [osmN.x[i] for i in traj]
        poiy = [osmN.y[i] for i in traj]
        ax.plot(poix,poiy, c='red',marker='o', linestyle='dashed')
        for an in x1: 
            ax.annotate(an,poi[an])
        
    ax.plot(lon,lat,c='green',marker='*', linestyle=' ')
    x2 = np.arange(0,len(pointsp))
    for an in x2: 
        ax.annotate(an,pointsp[an])
    plt.show()

def PlotOSMids(traj,ogp):
    
    x1 = np.arange(0,len(traj))
    osmN,osmE = ox.graph_to_gdfs(ogp,edges=True,nodes=True)
    poi = [(osmN.x[i],osmN.y[i]) for i in traj]
    poix = [osmN.x[i] for i in traj]
    poiy = [osmN.y[i] for i in traj]
    
    fig, ax = ox.plot_graph(ogp, show=False, close=False)
    ax.plot(poix,poiy, c='red',marker='o', linestyle='dashed')
    for an in x1: 
        ax.annotate(an,poi[an])
    plt.show()

def CreateBbox(listLon,listLat):
    lon = listLon
    lat = listLat
    #points = list(zip(lon,lat))
    bbox = [max(lat),min(lat), max(lon), min(lon)] ##[maxLat,minLat,maxLon,minLon]
    og = ox.graph_from_bbox(bbox[0], bbox[1], bbox[2], bbox[3],truncate_by_edge=True,custom_filter=drive_filter)
    ogp = ox.projection.project_graph(og)
    return ogp

def osmTodis(o1,o2,og):#,omsN,osmE): # takes two osmids, return two lists, distance time velocity etc.
    #tl = [t1,t2]
    #tdelta = datetime.strptime(t2, FMT) - datetime.strptime(t1, FMT)
    try:
        path = nx.shortest_path(og,o1,o2,weight='length')
    except:
        path = []
    if path:dist = sum([og[path[i]][path[i+1]][0]['length'] for i in range(len(path)-1)])
    else: dist = 1000000 ## almost infinite for no path
    #vel = []
    return path,dist


def OsmRouteFlatten(list_osm_points,ogp):
    #osmNF =
    #osmEF =
    OsmList =  [[]] * len(list_osm_points)
    for idx,item in enumerate(list_osm_points):
        path = []
        for idy,item1 in enumerate(item[1::]):
            o1 = item[idy]
            o2 = item1
            path = path + osmTodis(o1,o2,ogp)[0]
            #print(path)
        OsmList[idx] = path
    return OsmList


stime = []
fname = []
for tiCount in range(1,97):
    M = (tiCount * 15 )%60
    H = (tiCount * 15) //60
    timeM = '%02d'%M
    timeH = '%02d'%H
    timeStr = str(timeH)+ ':' + str(timeM) + ':00'
    stime.append(timeStr)
    

#ODM = pd.read_csv('ODFinal.csv',index_col=0)
#OD = pd.read_csv('ScatLinkDataSmall.csv',index_col=0)

AF = pd.read_csv('ODFinal.csv',index_col=0)
linksOD = [ast.literal_eval(item) for item in AF.index.to_list()]
AF = AF.loc[:, (AF != 0).any(axis=0)]
nameOD = [ast.literal_eval(item) for item in AF.dtypes.index.to_list()] #list(AF.dtypes.index[1::])

A = AF.values
NumColA = len(A[0])
NumRowA = len(A)


#B = np.zeros((NumRowA,NumColA))
LC = pd.read_csv('LinkCount_2019.csv',index_col=0)
#LC = LC.iloc[0:10472] ##UNCOMMENT THIS LINE FOR OTHER YEARS>>>> ONLY TRUE FOR 2016 SCATS DATA

LC['Link'] = list(zip(LC.LinkFrom.to_list(),LC.LinkTo.to_list()))
LC['Date'] = pd.to_datetime(LC['Date'])# ('Date',drop=True,inplace=True)
LCwk = LC[LC.Date.dt.dayofweek <5]
#start_date = '2016-01-01' 
#end_date = '2016-03-31'

Y_all_time = pd.DataFrame(columns=stime)
iar = pd.DataFrame(index=nameOD,columns=stime)
iar.fillna(0,inplace=True)
#iar.apply(lambda x: print(x.index),axis=1)
#iar.index.to_series().apply(lambda x: print(x))
#a = LCwk[LCwk.Link == x][stime].T.values.tolist() ## weekdays each 15 mins counts
#iar.index.to_series().apply(lambda x: LCwk[LCwk.Link == x][stime].T.values.tolist())
#iaR = iar.copy()
for idx in linksOD:
    temp = pd.DataFrame([LCwk[LCwk.Link == idx][stime].T.values.tolist()],columns=stime) #,index=idx)
    Y_all_time = Y_all_time.append(temp,ignore_index=True)
Y_all_time.index = linksOD
    #iaR[iar.index==idx]
    #LCwk[LCwk.Link == idx][stime].T.values.tolist()
    #iaR[iar.index==idx][stime] = pd.DataFrame(temp)
    #iaR[iar.index==idx][stime] = LCwk[LCwk.Link == idx][stime]
## create Y for each time interval in Y_all_time

Y_all_time.to_csv('Y_all_time_2019.csv')

for idx,tt in enumerate(stime): #each colum of Y_all iterate->
    Lam = (np.ones(NumColA))
    LamAll = np.ones(NumColA)
    #tempo = Y_all_time[tt]
    #Y = tempo.apply(lambda x: np.array(x).T).to_numpy()
    Y = np.asarray(Y_all_time[tt].values.tolist())  #Link.values[:,24::]
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

    AFinal = np.vstack((A,F_B))
    YFinal = np.vstack((Y_bar,F_S))

    #AFinal = A
    #YFinal = Y_bar

    TotalRow = len(YFinal)
    TempLam = np.ones(NumColA)
    #Y[0:4:3, 0:4:3]
    for numIter in range (0,50):
        for Countj in range (0,NumColA):
            #print(Countj)
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
    print(tt)
    iar[tt] = TempLam
    
#Lall = pd.DataFrame(LamAll)
name = '2019_IAR.csv' 
iar.to_csv(name)

