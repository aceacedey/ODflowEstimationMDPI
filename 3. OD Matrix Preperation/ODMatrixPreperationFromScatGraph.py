import math
import pandas as pd
import numpy as np
import networkx as nx
import re
import matplotlib.pyplot as plt
import scipy
from scipy.special import comb
import itertools
from numpy.linalg import matrix_power

from itertools import permutations

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
      #print(SumQ)
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
   print('----------------------- LAST COL-----------')
   ODFinal[ALI,i] = Res[0:-1]
   return ODFinal 
         

def StateTransitionMatrix(q0,q1,m):
   Q0 = q0
   Q1 = q1
   M = m
   


G = nx.read_edgelist("testEsmall.txt", create_using=nx.DiGraph(),data=(('weight',float),))

d = {}
with open("nodesXY.txt") as f:
    for line in f:
       key, *val=line.split()
       d[key] = [float(val[0]),float(val[1])]

nx.draw(G, pos=d, with_labels = True, node_size= 3000, node_color='r', alpha=.8,width=3 )
#A = nx.adjacency_matrix(G, nodelist=None, weight='weight')
#edge_labels = nx.get_edge_attributes(G,'weight')
#nx.draw_networkx_edge_labels(G, pos=d, labels = edge_labels)
#plt.Figure()
#plt.show()


NodeList = list(G.node())

#results = itertools.combinations(NodeList,2) # convert the combination iterator into a numpy array
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
      #print(Q11)
   else:
      print(Ori)
      print(Dest)
      

indexR = list(G.edges())

ODdataFrame = pd.DataFrame(ODFinal,index=indexR)
#ODdataFrame = ODdataFrame.loc[:, (ODdataFrame!= 0).any(axis=0)]
#headera = (list(itertools.compress(header, (ODdataFrame!= 0).any(axis=0))))
ODdataFrame.to_csv('ODFinal2.csv',index=index,header=headerOD)
