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
from itertools import permutations

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



def CalProba(SNode,DNode):
   U = SNode
   V = DNode
   if (U!=V):
      SP = list(nx.all_shortest_paths(G, source=U, target=V)) 
      SPLen = len(SP)
      nSP = np.array(SP)
      
      print('Shortest Paths' + str(SP))
      #for CountSP in range (0,SPLen):
      print('----' + str(SPLen) + ' ' + str(U) + ' ' + str(V) + '----')
      AllSuc = list(G.successors(U))
      CurrShortestPathSucc = list(nSP[:,1])
      #print('Cur Shortest Path ' + str(CurrShortestPathSucc))
      TempSuc = list(set(CurrShortestPathSucc).intersection(set(AllSuc)))
      #print('Intersection with Succ and Shortest Path: ' + str(TempSuc))
      NumTempSuc = len(TempSuc)   
      prob = 1/NumTempSuc
         
         #print('Prob: ' + str(prob))
      for CountTempSuc in range (0,NumTempSuc):
         
         TempDest = TempSuc[CountTempSuc]
        
         CurEdge = np.array(list((U,TempDest)))
         c1 = np.where(CurEdge[0] == ODRow[:,0])
         c2 = np.where(CurEdge[1] == ODRow[:,1])
         indexRow = np.intersect1d(c1,c2)
         OD[indexRow,i] = prob
         #Q_0[indexRow,State] = prob
         if TempDest == V:
            return 
         else:
            #State = State + 1
            CalProba(TempDest,V)

   else:
      return


def BuildQ1(So,De,SL,Max,Idx):
   StateList = SL
   #print('Active Link Index:')
   #print(StateList)
   M = Max
   #print(M)
   S = So
   D = De
   ALI = np.array(Idx)
   ALI = np.array(Idx)
   Q0 = np.zeros(M)
   Q0index = np.where(StateList[:,0]==S)
   Q0[Q0index] = OD[ALI[0,Q0index],i]
   print('Q0 ---->')
   #print(Q0)
   for count in range(0,M-1):
      PrevStateSource = StateList[count,0]
      PrevStateSink = StateList[count,1]
      if PrevStateSink == D:
         Q1[count,-1] = 1
         #print('here..............')
      else:
         for count2 in range(0,M-1):
            NextStateSource = StateList[count2,0]
            NextStateSink = StateList[count2,1]
            if PrevStateSink == NextStateSource:
               #print('here2..............')
               Q1[count,count2] = OD[ALI[0,count2],i]
   #print(Q1)
   #StateTransitionMatrix(Q0,Q1,M)C
   Q1[-1,-1] = 1
   MulCount = 0
   Q_temp = Q1
   #print(Q_temp)
   SumQ = np.zeros((M,M))
   Res = np.zeros(M)
   numChain = 0
   #for numChain in range (0,M):
   while np.count_nonzero(SumQ[:,-1]) != M: 
      SumQ = matrix_power(Q_temp, numChain)
      Res = Res + np.matmul(Q0,SumQ)
      numChain = numChain + 1
     # print(SumQ)
   #Res = np.matmul(Q0,SumQ)
   #print(Q_temp[:,-1])
   #print(np.count_nonzero(Q_temp[:,-1]))
   #while np.count_nonzero(Q_temp[:,-1]) != M: 
      #Res = Res + np.matmul(Res,Q_temp) 
      #Q_temp = np.matmul(Q_temp,Q_temp)
      #print(i)
   #print('----------------------- LAST COL-----------')
   #print(SumQ[:,-1])
   #print(Res[0:-1])
   #print('----------------------- LAST COL-----------')
   ODFinal[ALI,i] = Res[0:-1]
   return ODFinal 
         

def StateTransitionMatrix(q0,q1,m):
   Q0 = q0
   Q1 = q1
   M = m

inProj = Proj(init='epsg:4326')     # original coordinate system
outProj = Proj(init='epsg:28355') 

#slink = pd.read_csv("ScatLinkDataSmall.csv",index_col=0)
ut1 = pd.read_csv('melbourne-tz-2019-1-OnlyWeekdays-HourlyAggregate.csv',index_col=0)
ut1.reset_index(inplace=True) ## important
filename = 'melbourne-tz-2019-1-OnlyWeekdays-HourlyAggregate'
scat = pd.read_csv('ScatsLocations_100m.csv',index_col=0)

scatXY = pd.read_csv('scatXY.csv')
position_dict = scatXY.set_index('Scatid').T.to_dict('list')
scatFT = pd.read_csv('scatFromTo.csv')
G = nx.from_pandas_edgelist(scatFT, 'From', 'To','w',create_using=nx.DiGraph())
scatXY = scatXY[scatXY.Scatid.isin(list(G.nodes))]

### Drawing ScatGraph here
##nx.draw(G, pos=position_dict, with_labels = True, node_size= 700, node_color='white', alpha=.8,width=1,font_size=14)
#A = nx.adjacency_matrix(G, nodelist=None, weight='weight')
#edge_labels = nx.get_edge_attributes(G,'weight')
#nx.draw_networkx_edge_labels(G, pos=position_dict, labels = edge_labels)
##plt.Figure()
##plt.show()

startTime = datetime.datetime.now()
NodeList = list(G.node())

#results = itertools.combinations(NodeList,2) # convert the combination iterator into a numpy array
permutes = itertools.permutations(NodeList,2)

ODCol = list(permutes)
SizeODCol = len(ODCol)
ODRow = list(G.edges())
SizeODRow = len(ODRow)
OD = np.zeros((SizeODRow,SizeODCol))
#Q11 = np.zeros((SizeODRow,SizeODCol))
ODFinal = np.zeros((SizeODRow,SizeODCol))
colNames = ['O', 'D',  'hod', 'ttmean', 'ttstd', 'geomean', 'geostd']

scatTT = pd.DataFrame(columns = colNames)
#scatRoute = pd.DataFrame(columns = ['O','D','ScatRoute','OsmRoute'])


if 1 == 1:                                 
    inProj = Proj(init='epsg:4326')     # original coordinate system
    outProj = Proj(init='epsg:28355')  
    lonp, latp = transform(inProj, outProj, list(scatXY.lon.values), list(scatXY.lat.values))
    pointsp = list(zip(lonp,latp))
    #bbox = [-37.7962, -37.8261, 144.9904, 144.9337]
    lon, lat = list(scatXY.lon.values), list(scatXY.lat.values)
    bbox = [max(lat),min(lat), max(lon), min(lon)]
    drive_filter = ox.core.get_osm_filter('drive') 
    #og = ox.graph_from_bbox(bbox[0], bbox[1], bbox[2], bbox[3],truncate_by_edge=True,custom_filter=drive_filter,simplify=True)
    #og = ox.graph_from_bbox(bbox[0], bbox[1], bbox[2], bbox[3],truncate_by_edge=True,bidirectional=True,retain_all=True,simplify =False)
    og = ox.graph_from_file('Melbourne.osm',bidirectional=False,retain_all=False,simplify =False)

    #fullogp = ox.projection.project_graph(fullog, to_crs=28355)
    ogp = ox.projection.project_graph(og, to_crs=28355) ##
    osmN = ox.graph_to_gdfs(ogp,edges=False)
    osmE = ox.graph_to_gdfs(ogp,nodes=False)
    AllScatPts = MultiPoint(osmN.geometry.to_list())
    
    fig, ax = ox.plot_graph(ogp, show=False, close=False)
    plt.show()



for item in ODCol:
    ori = item[0]
    des = item[1]
    oriU = scat[scat.SITE_NO == int(ori)].Uber_id
    oriU = ast.literal_eval(oriU.to_list()[0]) ## list of uber movement zones as Origin for one scat id
    oriU = list(np.array(oriU,dtype=int))
    desU = scat[scat.SITE_NO == int(des)].Uber_id
    desU = ast.literal_eval(desU.to_list()[0]) ## list of uber movement zones as Destination for one scat id
    desU = list(np.array(desU,dtype=int))

    ##ut1[ut1.sourceid.isin(oriU)]
    ## for each source id calculate destination average mean and average standard deviation
    tempOD = ut1[ut1.sourceid.isin(oriU) & ut1.dstid.isin(desU)]
    ## groupby tempOD for each hour of a day

    tempgrp = tempOD.groupby(['hod']).groups
    #length_of_subrows = len(tempgrp.keys())
    #Pvalues = pd.DataFrame(length_of_subrows * [[ori,des,filename]])
    ## for each hour, calculate mean and standard deviation ....
    for tk in tempgrp: 
        loi = tempgrp[tk] 
        #print(loi)
        
        ttMean = tempOD.loc[loi,'mean_travel_time'].mean()
        ttStd = tempOD.loc[loi,'standard_deviation_travel_time'].mean()
    
        gmean = tempOD.loc[loi,'geometric_mean_travel_time'].mean()
        gstd = tempOD.loc[loi,'geometric_standard_deviation_travel_time'].mean()
        # create a dataframe ['ScatSource', 'ScatDest', 'hod', 'ttmean', 'ttstd', 'description']
        scatTT = scatTT.append(pd.DataFrame([[ori,des,tk,ttMean,ttStd,gmean,gstd]],columns=colNames),ignore_index=True)
        
scatTT.to_csv('scatTT_Q1_2019.csv')



######## comment/Uncomment the following if ODFinal is not found in right format or Unsure # dated 18/05/2020

##
##if 1==1: ## Plot SCAT ids along with nearest osmids   
##    fig, ax = ox.plot_graph(ogp, show=False, close=False)
##    ax.plot(lonp,latp,c='green',marker='o', linestyle=' ')
##    for i in range(len(scatXY)):
##        ax.annotate(int(scatXY.iloc[i].Scatid),pointsp[i])
##    x = []
##    y = []
##    pts = MultiPoint(osmN.geometry.to_list())
##    for item in list_scat_points[0]:
##        
##        #n1 = ox.get_nearest_node(ogp,item)
##        #n1 = ox.get_nearest_node(ogp,(item[1],item[0]))
##        
##        n = nearest_points(Point(item),pts)
##        n1 = osmN[osmN.geometry == n[1]].osmid.values[0]
##        x.append(osmN.x[n1])
##        y.append(osmN.y[n1])
##    #n2 = ox.get_nearest_node(ogp,list_scat_points[0][-1])
##    ax.plot(x,y,c='red',marker='o', linestyle='-')
##
##    
##    x1 = np.arange(0,len(OsmRoute[0]))
##    traj = OsmRoute[0]
##    poi = [(osmN.x[i],osmN.y[i]) for i in traj]
##    poix = [osmN.x[i] for i in traj]
##    poiy = [osmN.y[i] for i in traj]
##    
##    #fig, ax = ox.plot_graph(ogp, show=False, close=False)
##    ax.plot(poix,poiy, c='blue',marker='o', linestyle='dashed')
##    for an in x1: 
##        ax.annotate(an,poi[an])
##    #ax.plot(osmN.x[n2],osmN.y[n2],c='red',marker='o', linestyle=' ')
##    plt.show()
##    
##if 1 == 1:
##    OsmList =  [[]] * len(list_osm_points)
##    for idx,item in enumerate(list_osm_points):
##        path = []
##        for idy,item1 in enumerate(item[1::]):
##            o1 = item[idy]
##            o2 = item1
##            path = path + osmTodis(o1,o2,ogp)[0]
##            print(path)
##        OsmList[idx] = path
##
##NodeList = list(G.node())
##
##ODCol = list(permutations(NodeList,2))
##headerOD = ODCol
##SizeODCol = len(ODCol)
##
##ODRow = np.array(list(G.edges()))
##SizeODRow = len(list(G.edges()))
##OD = np.zeros((SizeODRow,SizeODCol))
##Q11 = np.zeros((SizeODRow,SizeODCol))
##ODFinal = np.zeros((SizeODRow,SizeODCol))
##
##for i in range(0,SizeODCol):
##   Ori = ODCol[i][0]
##   Dest = ODCol[i][1]
##   print(i)
##   print(Ori)
##   print(Dest)
##   #print(Dest)
##   State = 0
##   
##   if nx.has_path(G, Ori, Dest):
##      #iSP = list(nx.all_shortest_paths(G, source=U, target=V)) 
##      #iSPLen = len(iSP)      
##      #Q1 = np.zeros   
##      CalProba(Ori,Dest)
##      #print(OD[:,i])
##      activeLinkIndex = np.nonzero(OD[:,i])
##      MaxState = np.array(activeLinkIndex).size + 1
##      Q0index = np.nonzero(OD[:,i])
##      Q1 = np.zeros((MaxState,MaxState))
##      #Q11 = np.zeros((MaxState))
##      #Q0 = 
##      QCol = np.array(ODRow[activeLinkIndex])  ### Pass edges with positive probability
##      #print('try......')
##      #print(Q0index)
##      Q11 = BuildQ1(Ori,Dest,QCol,MaxState,activeLinkIndex)
##      
##      print(Q11)
##   else:
##      print(Ori)
##      print(Dest)
##   
##indexR = list(G.edges())
##ODdataFrame = pd.DataFrame(ODFinal,index=indexR)
##ODdataFrame.to_csv('ODFinal.csv',index=index,header=headerOD)

#a = pd.read_csv('ODFinal.csv',index_col=0)
#[ast.literal_eval(item) for item in a.index.to_list()]
## path = [3208518282, 955611524, 1985546057, 137182984, 471341019, 477075428, 6775928319, 3933090382, 30943325, 1985551132, 1985549834]



##    U = ori
##    V = des
##    #list_scat_points = []
##    if (U!=V) and nx.has_path(G,U,V) == True:
##        SP = list(nx.all_shortest_paths(G, source=U, target=V)) 
##        SPLen = len(SP)
##        list_scat_points = np.zeros(np.shape(SP),dtype=object) #[[]] * SPLen
##        list_osm_points = np.zeros(np.shape(SP),dtype=object)
##        nSP = np.array(SP)
##        for idx,tempSP in enumerate(SP):
##
##            for idy,tempNode in enumerate(tempSP):
##                temp = scatXY[scatXY.Scatid == tempNode]
##                xp,yp = transform(inProj, outProj,temp.lon.values[0] , temp.lat.values[0])
##                
##                list_scat_points[idx][idy] = (xp,yp)
##                n = nearest_points(Point(xp,yp),AllScatPts)
##                temp_osmid = osmN[osmN.geometry == n[1]].osmid.values[0]
##                #ox = osmN.x[temp_osmid])
##                #oyy.append(osmN.y[temp_osmid])
##                list_osm_points[idx][idy] = temp_osmid
##
##    OsmRoute = OsmRouteFlatten(list_osm_points,ogp) ### find intermediate osmids
##    #PlotOSMids(OsmRoute[0],ogp)
##    #break
##    tempdf = pd.DataFrame([[ori,des,SP,list_scat_points,OsmRoute]],columns = ['O','D','ScatRoute','PointRoute','OsmRoute'])
##    scatRoute = scatRoute.append(tempdf,ignore_index=True)
## 
##    
##scatRoute.to_csv('ScatOsmPointRoutes.csv')

