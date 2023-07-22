import shapely
from shapely.geometry import shape, Point, mapping
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



def PlotTraj(traj,pointsp,ogp):
    lon = []
    lat = []
    for item in pointsp:
        lon.append(item[0])
        lat.append(item[1])
        
    x1 = np.arange(0,len(traj))
    x2 = np.arange(0,len(pointsp))
    osmN,osmE = ox.graph_to_gdfs(ogp,edges=True,nodes=True)
    poi = [(osmN.x[i],osmN.y[i]) for i in traj]
    poix = [osmN.x[i] for i in traj]
    poiy = [osmN.y[i] for i in traj]
    
    fig, ax = ox.plot_graph(ogp, show=False, close=False)
    ax.plot(poix,poiy, c='red',marker='o', linestyle='dashed')
    for an in x1: 
        ax.annotate(an,poi[an])
    
    ax.plot(lon,lat,c='green',marker='*', linestyle=' ')
    for an in x2: 
        ax.annotate(an,pointsp[an])
    plt.show()


def CreateBbox(listLon,listLat):
    lon = listLon
    lat = listLat
    #points = list(zip(lon,lat))
    bbox = [max(lat),min(lat), max(lon), min(lon)] ##[maxLat,minLat,maxLon,minLon]
    og = ox.graph_from_bbox(bbox[0], bbox[1], bbox[2], bbox[3],truncate_by_edge=True,custom_filter=drive_filter)
    ogp = ox.projection.project_graph(og)
    return ogp

#def SCATAG(scatid,df):# enter a scat id, get aggregated counts 
G = nx.read_edgelist("testEsmall.txt", create_using=nx.DiGraph(),data=(('weight',float),))
slink = pd.read_csv("ScatLinkDataSmall.csv",index_col=0)
#uz = pd.read_csv('Melbourne_movement_zones_Uber.csv',index_col=0)
scat = pd.read_csv('ScatsLocations_100m.csv',index_col=0)


d = {}
with open("nodesXY.txt") as f:
    for line in f:
       key, *val=line.split()
       d[key] = [float(val[0]),float(val[1])]
nx.draw(G, pos=d, with_labels = True, node_size= 700, node_color='white', alpha=.8,width=1,font_size=14)
#A = nx.adjacency_matrix(G, nodelist=None, weight='weight')
edge_labels = nx.get_edge_attributes(G,'weight')
#nx.draw_networkx_edge_labels(G, pos=d, labels = edge_labels)
plt.Figure()
plt.show()

startTime = datetime.datetime.now()

inProj = Proj(init='epsg:4326')     # original coordinate system
outProj = Proj(init='epsg:28355')  


bbox = [-37.7962, -37.8261, 144.9904, 144.9337]
drive_filter = ox.core.get_osm_filter('drive') 
og = ox.graph_from_bbox(bbox[0], bbox[1], bbox[2], bbox[3],truncate_by_edge=True,custom_filter=drive_filter,simplify=True)

#fullog = ox.graph_from_bbox(bbox[0], bbox[1], bbox[2], bbox[3],truncate_by_edge=True)

ogp = ox.projection.project_graph(og, to_crs=28355) ##
osmN = ox.graph_to_gdfs(ogp,edges=False)
osmE = ox.graph_to_gdfs(ogp,nodes=False)
fig, ax = ox.plot_graph(ogp, show=False, close=False)
plt.show()



##
##
##
##
##buf = 100 # in meter
##
##scat['BBox'] = scat.apply(lambda x: box(x.X - buf, x.Y - buf, x.X + buf, x.Y + buf),axis=1)
##
##uz['geometryP'] = uz['geometryP'].apply(wkt.loads) ## conevert string to polygon back
##
##uzGeo = gpd.GeoSeries(uz.geometryP)  ## create a geoseries
##
###uz.MOVEMENT_ID[a.intersects(p2)].to_list()
##
##scat['Uber_id'] = scat.BBox.map(lambda x: uz.MOVEMENT_ID[uzGeo.intersects(x)].to_list())
##
##
##
##
##
##### Example polygon
##xy = [[130.21001, 27.200001], [129.52, 27.34], [129.45, 27.1], [130.13, 26.950001]]
##polygon_shape = Polygon(xy)
##
##
##p = box(129.5, -27.0, 129.75, 27.25)  #  (minx, miny, maxx, maxy) or
##x = 100
##y= 20 
##p1 = box(x - buf,y - buf,x + buf,y + buf)
##### The intersection
##print(polygon_shape.intersection(p1).area)
### Example grid cell
##if 1 == 1:
##    plt.plot(scat.X.values[0],scat.Y.values[0],'bo')
##    p1 = scat.BBox.values[0]
##    plt.plot(*p1.exterior.xy)
##    plt.show()
#p2 = scat[scat.SITE_NO == 2023].BBox.values[0]    
#p2.intersects(uz.geometryP.any())


#poly.geometry.map(lambda x: x.intersects(line.geometry.any()))
#uz.geometryP.apply(shapely.wkt.loads)
## for each point in scat, scan through 
##
##
##
##
##from shapely.ops import cascaded_union
##from rtree import index
##idx = index.Index()
##
### Populate R-tree index with bounds of grid cells
##for pos, cell in enumerate(grid_cells):
##    # assuming cell is a shapely object
##    idx.insert(pos, cell.bounds)
##
### Loop through each Shapely polygon
##for poly in polygons:
##    # Merge cells that have overlapping bounding boxes
##    merged_cells = cascaded_union([grid_cells[pos] for pos in idx.intersection(poly.bounds)])
##    # Now do actual intersection
##    print(poly.intersection(merged_cells).area)
## 
###startTime = datetime.now()
##FMT = '%H:%M:%S'
##
