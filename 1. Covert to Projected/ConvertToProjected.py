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

#def SCATAG(scatid,df):# enter a scat id, get aggregated counts 
    
f = open('melbourne_tz.json') 
  
# returns JSON object as  
# a dictionary 
data = json.load(f) 
Uzones = json_normalize(data['features'])

crs = {'init': 'epsg:4326'}

gdf = gpd.GeoDataFrame.from_features(data["features"], crs={'init': 'epsg:4326'})
gdfp = gdf.copy()

#>>> gdf = gpd.GeoDataFrame.from_features(data["features"], crs={'init': 'epsg:4326'})
		     
#>>> gdfp = gdf.copy()
		     
gdf_proj = gdfp['geometry'].to_crs(epsg=28355)
gdfp['geometryP'] = gdf_proj
gdfp.to_csv('Melbourne_movement_zones_Uber.csv')

#Uzones['geometry.coordinates'].apply(lambda x: x[0],axis=1)

scat = pd.read_csv('AllScatIDLocations.csv')
pointsp = scat.apply(lambda x: utm.from_latlon(x.POINT_Y,x.POINT_X)[0:2],axis=1) ## utm package take lat,lon and convert to lon, lat
scat['X'] = scat.apply(lambda x: utm.from_latlon(x.POINT_Y,x.POINT_X)[0],axis=1)
scat['Y'] = scat.apply(lambda x: utm.from_latlon(x.POINT_Y,x.POINT_X)[1],axis=1)

from shapely.geometry import box, Polygon

buf = 50 # in meter

scat['BBox'] = scat.apply(lambda x: box(x.X - buf, x.Y - buf, x.X + buf, x.Y + buf),axis=1)

#uz['geometryP'] = uz['geometryP'].apply(wkt.loads) ## conevert string to polygon back
uz = gdfp
uzGeo = gpd.GeoSeries(uz.geometryP)  ## create a geoseries

#uz.MOVEMENT_ID[a.intersects(p2)].to_list()

scat['Uber_id'] = scat.BBox.map(lambda x: uz.MOVEMENT_ID[uzGeo.intersects(x)].to_list())

scat.to_csv('ScatsLocations.csv')

##
##from shapely.geometry import box, Polygon
### Example polygon 
##xy = [[130.21001, 27.200001], [129.52, 27.34], [129.45, 27.1], [130.13, 26.950001]]
##polygon_shape = Polygon(xy)
### Example grid cell
##grid_cells = box(129.5, -27.0, 129.75, 27.25)
### The intersection
##polygon_shape.intersection(grid_cells).area
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


##filesa2 = "1270055001_sa2_2016_aust_shape\SA2_2016_AUST.shp"
##filesa3 = "1270055001_sa3_2016_aust_shape\SA3_2016_AUST.shp"
###filename = "SA2_2016_AUST.shape"
###shape = fiona.open(file)
##
##import geopandas as gpd
##
##sa2 = gpd.read_file(filesa2)
##sa3 = gpd.read_file(filesa3)
##
##scat = pd.read_csv("AllScatIDLocations.csv")
##Melb_sa3 = sa3[sa3.SA3_NAME16=='Melbourne City'].geometry.values[0]
#### for each unique SA3_MAIN16, do the sameSA3 
##melsa2 = sa2[sa2.SA3_NAME16=='Melbourne City'] ## this dataframe has all sa2 polygons under SA3_melbourne
##Melbourne_SA3_code = melsa2.SA3_CODE16.values[0] #melbourne city sa3 code 16
###point = shapely.geometry.Point(144.962178,-37.812171) ## Longitude first, a random point in CBD
##
###nameofSA2 = melsa2.SA2_NAME16[melsa2.geometry.contains(point)].values[0] #melsa2[melsa2.geometry.contains(point)] ## gives the name of the SA2 where this point contains 
##ListSA3 = melsa2.SA2_MAIN16.to_list() ## List of SA3 zones under nameofSA2 
##
##tp = shapely.geometry.Point ## object to create a point
##
##scat['SA2'] = scat.apply(lambda x: melsa2.SA2_MAIN16[melsa2.geometry.contains(tp(x.POINT_X,x.POINT_Y))].values,axis=1)
###scat.SA2.values.size == 0
##
##
##scatdf = scat[scat.SA2.apply(lambda x: x.size != 0)] ### a dataframe with scatid and their correspoinding sa2
##newSA2 = scatdf.SA2.apply(lambda x: x[0])
##scatsa2df = scatdf.replace(scatdf.SA2.values,newSA2.values)
##
##mypath = "C:\Work Temporary Storage\SCATSRawData\VSDATA_ALL"
##
##dataframes = ''
##filenames = [f for f in listdir(mypath) if isfile(join(mypath, f))]
##numFiles = len(filenames)
###Rows = len(SLDS)
##timeintervals = 96
###YhatFinal = np.zeros((timeintervals,numFiles,Rows))
##Header = []
##i = 0
##date_list = []
##sa3_scat = dict()#{'name': [],'scat_id':[]}
##
##for item in ListSA3:
##    temp = scatsa2df[scatsa2df.SA2 == item]
##    print(item)
##    #value = temp.SITE_NO.to_list()
##    sa3_scat.update({item:temp.SITE_NO.to_list()})
##
##sa3df = pd.DataFrame(index=range(numFiles),columns=ListSA3) ## wide format
##sa3df.fillna(0,inplace=True)
##
##### lets create a dataframe in long format
##total_rows = numFiles * len(ListSA3)
##
##
##vtime = [] #['V00', 'V01', 'V02', 'V03', 'V04', 'V05', 'V06', 'V07', 'V08', 'V09', 'V10', 'V11', 'V12', 'V13', 'V14', 'V15', 'V16', 'V17', 'V18', 'V19', 'V20', 'V21', 'V22', 'V23', 'V24', 'V25', 'V26', 'V27', 'V28', 'V29', 'V30', 'V31', 'V32', 'V33', 'V34', 'V35', 'V36', 'V37', 'V38', 'V39', 'V40', 'V41', 'V42', 'V43', 'V44', 'V45', 'V46', 'V47', 'V48', 'V49', 'V50', 'V51', 'V52', 'V53', 'V54', 'V55', 'V56', 'V57', 'V58', 'V59', 'V60', 'V61', 'V62', 'V63', 'V64', 'V65', 'V66', 'V67', 'V68', 'V69', 'V70', 'V71', 'V72', 'V73', 'V74', 'V75', 'V76', 'V77', 'V78', 'V79', 'V80',
###'V81', 'V82', 'V83', 'V84', 'V85', 'V86', 'V87', 'V88', 'V89', 'V90', 'V91', 'V92', 'V93', 'V94', 'V95']
##stime = []
##for tiCount in range(1,97):
##    num = '%02d'%(tiCount-1)
##    Vs = 'V'+str(num)
##    vtime.append(Vs)
##    M = (tiCount * 15 )%60
##    H = (tiCount * 15) //60
##    timeM = '%02d'%M
##    timeH = '%02d'%H
##    timeStr = str(timeH)+ ':' + str(timeM) + ':00'
##    stime.append(timeStr)
##    
##wide_to_longCols = ['Date','SA3','SA2','VOLUME_24HOUR'] + stime
##
##sa3dfL = pd.DataFrame(index=range(total_rows),columns=wide_to_longCols ) 
##sa3dfL.fillna(0,inplace=True)
##total_row_index = 0
##y0 = '2018' ## starting year of the VSDATA
##for idx,f in enumerate(filenames): ########
##    print(f)
##    list_vol = []
##    
##    y = f[7:11]
##    m = f[11:13]
##    d = f[13:15]
##    date = str(d) + '/' + str(m) + '/' + str(y)
##    date_list.append(date)
##    #print(date_list)
##    vsdf = pd.read_csv(join(mypath,f))
##    timeCols = vsdf.columns.to_list()[3:-4]
##    for k in sa3_scat:
##        #print(sa3_scat[k])
##        scat_to_search = sa3_scat[k]
##        
##        trimdf = vsdf[vsdf.NB_SCATS_SITE.isin(scat_to_search)] #vsdf[vsdf['NB_SCATS_SITE'] == 105] ##
##
##        sa3df.loc[idx,k] = trimdf.QT_VOLUME_24HOUR.sum() ### Add all counts for 24 hours 
##        # melbourne city sa3 number
##        sa3dfL.loc[total_row_index,'SA2'] = k
##        sa3dfL.loc[total_row_index,'Date'] = date
##        sa3dfL.loc[total_row_index,'VOLUME_24HOUR'] =  trimdf.QT_VOLUME_24HOUR.sum()
##        sa3dfL.loc[total_row_index,stime] =  trimdf[trimdf[vtime] >= 0][vtime].sum().values #trimdf[vtime].sum()
##        total_row_index =total_row_index + 1
##
##sa3dfL.loc[:,'SA3'] = [Melbourne_SA3_code] * total_rows
###sa3df.index = pd.to_datetime(date_list,format='%d/%m/%Y')
##sa3df.index = pd.to_datetime(date_list)#,format='%d/%m/%Y')
##sa3df.to_csv('b.csv')
##sa3dfL.to_csv('c.csv')
##
##
##
