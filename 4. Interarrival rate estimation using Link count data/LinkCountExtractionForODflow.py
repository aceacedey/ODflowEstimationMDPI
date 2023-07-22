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

import datetime
import ast

FMT = '%H:%M:%S'

import geopandas as gpd



SLDS = pd.read_csv('ScatLinkDataSmall.csv')
AllFrom = np.array(SLDS.From)
AllTo = np.array(SLDS.TO)
scat2SearchList = SLDS.SCATtoSearch.to_list()

## for each row of SLDS, find relevant lanes for a link

ODinitial = pd.read_csv('ODFinal.csv',index_col=0)
linksOD = [ast.literal_eval(item) for item in ODinitial.index.to_list()]

mypath = "C:\Work Temporary Storage\SCATSRawData\VSDATA_2016"

dataframes = ''
filenames = [f for f in listdir(mypath) if isfile(join(mypath, f))]
numFiles = len(filenames)
#Rows = len(SLDS)
timeintervals = 96
#YhatFinal = np.zeros((timeintervals,numFiles,Rows))
vtime = []
stime = []

for tiCount in range(1,97):
    num = '%02d'%(tiCount-1)
    Vs = 'V'+str(num)
    vtime.append(Vs)
    M = (tiCount * 15 )%60
    H = (tiCount * 15) //60
    timeM = '%02d'%M
    timeH = '%02d'%H
    timeStr = str(timeH)+ ':' + str(timeM) + ':00'
    stime.append(timeStr)
    
wide_to_longCols = ['Date','LinkFrom','LinkTo'] + stime

total_rows = len(SLDS) * len(filenames)#number of links * number of files
dfL = pd.DataFrame(index=range(total_rows),columns=wide_to_longCols ) 
dfL.fillna(0,inplace=True)
total_row_index = 0
Rows = len(SLDS)

for idx,f in enumerate(filenames): ########
    print(f)
    list_vol = []
    
    y = f[7:11]
    m = f[11:13]
    d = f[13:15]
    date = str(d) + '/' + str(m) + '/' + str(y)
    #date_list.append(date)
    #print(date_list)
    vsdf = pd.read_csv(join(mypath,f))
    timeCols = vsdf.columns.to_list()[3:-4]
    
    for linkCount in range(0,Rows):
        indexLinkAll = linkCount
        LinkToSearch = int(SLDS.SCATtoSearch[indexLinkAll])
        SensorNumbersL = list((SLDS.Lane1[indexLinkAll],SLDS.Lane2[indexLinkAll],SLDS.Lane3[indexLinkAll],SLDS.Lane4[indexLinkAll]))
        ActiveSensor = list(filter(None,SensorNumbersL))
        
        scat_to_search = LinkToSearch
        trimdf = vsdf[vsdf.NB_SCATS_SITE == scat_to_search] #trimdf = vsdf[vsdf.NB_SCATS_SITE.isin = scat_to_search] #vsdf[vsdf['NB_SCATS_SITE'] == 105] ##
        
        # melbourne city sa3 number
        dfL.loc[total_row_index,'LinkFrom'] = SLDS.From[linkCount] #str(SLDS.From[linkCount]) + '-' + str(SLDS.TO[linkCount])#[(SLDS.From[linkCount],SLDS.TO[linkCount])]
        dfL.loc[total_row_index,'LinkTo'] = SLDS.TO[linkCount]
        dfL.loc[total_row_index,'Date'] = date

        activedf = trimdf[trimdf.NB_DETECTOR.astype('int32').isin(ActiveSensor)] ### change column type
        dfL.loc[total_row_index,stime] = activedf[activedf[vtime] >= 0][vtime].sum().values # Find relevant lanes for that link and add those vertically
        total_row_index =total_row_index + 1
        
dfL['Link'] = list(zip(dfL.LinkFrom.to_list(),dfL.LinkTo.to_list()))
#sa3df.index = pd.to_datetime(date_list)#,format='%d/%m/%Y')
#fnamec = 'long_' + str(SA3Name) + '.csv'
year = "2016_Q1"
fnameb = 'LinkCount_' + str(year) + '.csv'
#sa3df.to_csv(fnamec)
dfL.to_csv(fnameb)
