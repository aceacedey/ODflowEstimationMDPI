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

plt.rcParams.update({'font.size': 24})

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
    if path:dist = [og[path[i]][path[i+1]][0]['length'] for i in range(len(path)-1)]
    else: dist = 1000000 ## almost infinite for no path
    #vel = []
    return path,dist


def OsmRouteFlatten(list_osm_points,ogp):
    #osmNF =
    #osmEF =
    #disList =  [[]] * len(list_osm_points)
    OsmList =  [[]] * len(list_osm_points)
    for idx,item in enumerate(list_osm_points):
        path = []
        dis = []
        for idy,item1 in enumerate(item[1::]):
            o1 = item[idy]
            o2 = item1
            path = path + osmTodis(o1,o2,ogp)[0]
            #dis.append(osmTodis(o1,o2,ogp)[1])
            #print(path)
        OsmList[idx] = path
        #disList[idx] = dis
    return OsmList#,disList

def DelShortEdge(pa,ogp): ## delete short edges in the beginging and end
    pa = pa
    dist = []
    #remove consecutive elements from path

##    for i in range(len(pa)-1):
##        d = osmE[(osmE.u == pa[i]) & (osmE.v == pa[i+1])]['length'].values[0]
    dist = [ogp[pa[i]][pa[i+1]][0]['length'] if pa[i] != pa[i+1] else 0 for i in range(len(pa)-1)]
            ####### filtering small edges in the osm route to reduce congestion: 
    for ele in dist:
        if ele <= 17: #17 meter as cuttoff distance 
            pa.remove(pa[0])
        else:
            break
    for ele in dist[::-1]:
        if ele <= 17:
            pa.remove(pa[-1])
        else:
            break
    return pa
def Scat2OsmDistance(P,ol):
    return [P.distance(osmN.loc[item].geometry) for item in ol] #, P is scat Point
def NearestOsmOfLink(ol,dlf,dlb):
    for ele in dlf:
        if ele <= 20: #17 meter as cuttoff distance 
            ol.remove(ol[0])
            dlb.remove(dlb[0])
        else:
            break
    for idx,ele2 in enumerate(dlb):
        if ele2 >= 20:
            continue
        else:
            break
    osmList = ol[0:idx]
    return osmList

inProj = Proj(init='epsg:4326')     # original coordinate system
outProj = Proj(init='epsg:28355') 

slink = pd.read_csv("ScatLinkDataSmall.csv",index_col=0)
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
nx.draw(G, pos=position_dict, with_labels = True, node_size= 1000, node_color='white', alpha=.8,width=1,font_size=20)
#A = nx.adjacency_matrix(G, nodelist=None, weight='weight')
#edge_labels = nx.get_edge_attributes(G,'weight')
#nx.draw_networkx_edge_labels(G, pos=position_dict, labels = edge_labels)
plt.Figure()
plt.show()

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
scatRoute = pd.DataFrame(columns = ['O','D','ScatRoute','OsmRoute'])


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
    AllOsmPts = MultiPoint(osmN.geometry.to_list())
    
    #fig, ax = ox.plot_graph(ogp, show=False, close=False)
    #plt.show()
scatXY['X'] = lonp
scatXY['Y'] = latp

scatXY['NPgeom'] = scatXY.apply(lambda x: nearest_points(Point(x.X,x.Y),AllOsmPts)[1],axis=1)#scatXY.apply(lambda x: x.
scatXY['NearestOsmid'] = scatXY.NPgeom.apply(lambda x: osmN[osmN.geometry == x].osmid.values[0])
if 1 == 1:
    scatFT['LinkID'] = scatFT.apply(lambda x: (int(x.From),int(x.To)),axis=1)
    scatFT['FromOsm'] = scatFT.From.apply(lambda x: scatXY.NearestOsmid[scatXY.Scatid == x].values[0])
    scatFT['ToOsm'] = scatFT.To.apply(lambda x: scatXY.NearestOsmid[scatXY.Scatid == x].values[0])
    scatFT['OsmPathLong'] = scatFT.apply(lambda x: osmTodis(x.FromOsm,x.ToOsm,ogp)[0],axis=1)
    #scatFT['OsmPath'] = scatFT.OsmPathLong.apply(lambda x: DelShortEdge(x,ogp))
    scatFT['FP'] = scatFT.From.apply(lambda x: Point(scatXY[scatXY.Scatid == x].X.values[0],scatXY[scatXY.Scatid == x].Y.values[0])) #,axis =1)
    scatFT['TP'] = scatFT.To.apply(lambda x: Point(scatXY[scatXY.Scatid == x].X.values[0],scatXY[scatXY.Scatid == x].Y.values[0]))
    scatFT['Forward'] = scatFT.apply(lambda x: Scat2OsmDistance(x.FP,x.OsmPathLong),axis=1)
    scatFT['Backward'] = scatFT.apply(lambda x: Scat2OsmDistance(x.TP,x.OsmPathLong),axis=1)
    #scatFT.apply(lambda x: np.array(x.OsmPathLong)[np.array(x.Forward)>10].tolist().index(np.array(x.OsmPathLong)[np.array(x.Forward)>10][0]) ,axis=1)
    scatFT['f'] = scatFT.apply(lambda x: np.array(x.OsmPathLong)[np.array(x.Forward) > 10].tolist(),axis=1)
    scatFT['b'] = scatFT.apply(lambda x: np.array(x.OsmPathLong)[np.array(x.Backward) < 10].tolist(),axis=1)
    scatFT['OsmPath'] = scatFT.apply(lambda x: x.f[0:x.f.index(x.b[0])+1],axis=1)
##scatFT.OsmPath.apply(lambda x: PlotOSMids(x,ogp))
scatFT.to_csv('scatFT.csv')
##    ol = scatFT.iloc[8].OsmPathLong
##    dlf =  scatFT.iloc[8].Forward
##    dlb = scatFT.iloc[8].Backward
##    f = np.array(ol)[np.array(dlf) > 10]
##    
##    b = np.array(ol)[np.array(dlb) < 10]
##    ix = f.tolist().index(b[0])
##    f.list()[0:ix]

for item in ODCol:
    ori = item[0]
    des = item[1]
    
    U = ori
    V = des
    #list_scat_points = []
    if (U!=V) and nx.has_path(G,U,V) == True:
        SP = list(nx.all_shortest_paths(G, source=U, target=V))
##        if len(SP) == 1:
##            sosm = scatXY.NearestOsmid[scatXY.Scatid == U].values[0]
##            eosm = scatXY.NearestOsmid[scatXY.Scatid == V].values[0]
##            path = nx.shortest_path(ogp,sosm,eosm,weight='length')
##            osmRoute = [DelShortEdge(path,ogp)] #### 
##            curScatPath = SP[0]
##            list_scat_points = [tuple(scatXY[scatXY.Scatid == item][['X','Y']].values[0]) for item in curScatPath]
##
##        else:
        list_osm_points = [[]] * len(SP) # np.zeros(np.shape(SP),dtype=object)
        list_links = [[]] * len(SP)
        for pathno,curPath in enumerate(SP):
            curPathdf = pd.DataFrame() ## create a dataframe to append relevant links in the scatPath 
            for idy,cNode in enumerate(curPath[1::]):
                link = (curPath[idy],cNode)  #links in curPath
                curPathdf = curPathdf.append(scatFT[scatFT.LinkID == link],ignore_index=True)
            list_osm_points[pathno] = curPathdf.OsmPath.sum()
            list_links[pathno] = curPathdf.LinkID.to_list()
        osmRoute = OsmRouteFlatten(list_osm_points,ogp)
##        for it in osmRoute:
##                PlotOSMids(it,ogp)
            
        tempdf = pd.DataFrame([[ori,des,SP,osmRoute]],columns = ['O','D','ScatRoute','OsmRoute'])
        scatRoute = scatRoute.append(tempdf,ignore_index=True)


scatRoute.to_csv('ScatOsmPointRoutes.csv')
#scatTT.to_csv('scatTT_Q1_2019.csv')
#scatRoute.Distance.apply(lambda x: [min(item) for item in x]).to_list()

######## comment/Uncomment the following if ODFinal is not found in right format or Unsure # dated 18/05/2020


if 1==1: ## Plot SCAT ids along with nearest osmids   
    fig, ax = ox.plot_graph(ogp, show=False, close=False)
    ax.plot(lonp,latp,c='green',marker='o', linestyle=' ')
    for i in range(len(scatXY)):
        ax.annotate(int(scatXY.iloc[i].Scatid),pointsp[i])
    x = []
    y = []
    pts = MultiPoint(osmN.geometry.to_list())
    for item in list_scat_points[0]:
        
        #n1 = ox.get_nearest_node(ogp,item)
        #n1 = ox.get_nearest_node(ogp,(item[1],item[0]))
        
        n = nearest_points(Point(item),pts)
        n1 = osmN[osmN.geometry == n[1]].osmid.values[0]
        x.append(osmN.x[n1])
        y.append(osmN.y[n1])
    #n2 = ox.get_nearest_node(ogp,list_scat_points[0][-1])
    ax.plot(x,y,c='red',marker='o', linestyle='-')

    
    x1 = np.arange(0,len(OsmRoute[0]))
    traj = OsmRoute[0]
    poi = [(osmN.x[i],osmN.y[i]) for i in traj]
    poix = [osmN.x[i] for i in traj]
    poiy = [osmN.y[i] for i in traj]
    
    #fig, ax = ox.plot_graph(ogp, show=False, close=False)
    ax.plot(poix,poiy, c='blue',marker='o', linestyle='dashed')
    for an in x1: 
        ax.annotate(an,poi[an])
    #ax.plot(osmN.x[n2],osmN.y[n2],c='red',marker='o', linestyle=' ')
    plt.show()
    
if 1 == 1:
    OsmList =  [[]] * len(list_osm_points)
    for idx,item in enumerate(list_osm_points):
        path = []
        for idy,item1 in enumerate(item[1::]):
            o1 = item[idy]
            o2 = item1
            path = path + osmTodis(o1,o2,ogp)[0]
            print(path)
        OsmList[idx] = path

NodeList = list(G.node())

ODCol = list(permutations(NodeList,2))
headerOD = ODCol
SizeODCol = len(ODCol)

ODRow = np.array(list(G.edges()))
SizeODRow = len(list(G.edges()))
OD = np.zeros((SizeODRow,SizeODCol))
Q11 = np.zeros((SizeODRow,SizeODCol))
ODFinal = np.zeros((SizeODRow,SizeODCol))

for i in range(0,SizeODCol):
   Ori = ODCol[i][0]
   Dest = ODCol[i][1]
   print(i)
   print(Ori)
   print(Dest)
   #print(Dest)
   State = 0
   
   if nx.has_path(G, Ori, Dest):
      #iSP = list(nx.all_shortest_paths(G, source=U, target=V)) 
      #iSPLen = len(iSP)      
      #Q1 = np.zeros   
      CalProba(Ori,Dest)
      #print(OD[:,i])
      activeLinkIndex = np.nonzero(OD[:,i])
      MaxState = np.array(activeLinkIndex).size + 1
      Q0index = np.nonzero(OD[:,i])
      Q1 = np.zeros((MaxState,MaxState))
      #Q11 = np.zeros((MaxState))
      #Q0 = 
      QCol = np.array(ODRow[activeLinkIndex])  ### Pass edges with positive probability
      #print('try......')
      #print(Q0index)
      Q11 = BuildQ1(Ori,Dest,QCol,MaxState,activeLinkIndex)
      
      print(Q11)
   else:
      print(Ori)
      print(Dest)
   
indexR = list(G.edges())
ODdataFrame = pd.DataFrame(ODFinal,index=indexR)
ODdataFrame.to_csv('ODFinal.csv',index=index,header=headerOD)

#a = pd.read_csv('ODFinal.csv',index_col=0)
#[ast.literal_eval(item) for item in a.index.to_list()]
## path = [3208518282, 955611524, 1985546057, 137182984, 471341019, 477075428, 6775928319, 3933090382, 30943325, 1985551132, 1985549834]

