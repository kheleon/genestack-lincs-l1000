import numpy as np
import pandas as pd
import os
import gzip
import csv

def getLevel2info():
    level2Files = []
    for file in os.listdir('./'):
        if file.endswith('gz') and 'Level2' in file:
            level2Files += ['./' + file]

    print ("reading Level2 files:",level2Files)

    i=0
    COI = {'id':0, 'det_plate':1, 'det_well':2, 'SM_LINCS_ID':3, 'SM_Pert_Type':4, 'SM_Time':5, 'SM_Time_Unit':6,
          'SM_Dose':7, 'SM_Dose_Unit':8}
    data = [[] for _ in xrange(len(COI))]
    indices = range(0,12)
    for filenum in range(len(level2Files)):
        with gzip.open(level2Files[filenum], 'r') as file:
            ind=0
            for line in file:
                tabline = line.rstrip().split('\t')
                if ind>20:
                    break
                ind=ind+1      
                if tabline[0] in COI.keys():
                    data[COI[tabline[0]]] = np.hstack([data[COI[tabline[0]]],tabline])

    exclude = [range(x,x+12) for x in np.where(data[0]=='id')[0]]
    exclude = [x for sublist in exclude for x in sublist]
    data = [data[i][np.setdiff1d(range(len(data[0])), exclude)] for i in range(0,len(data))]
    
    df = pd.DataFrame(data[0], columns=[COI.keys()[COI.values().index(0)]])
    for i in range(1, len(data)):
        df[COI.keys()[COI.values().index(i)]] = pd.Series(data[i], index=df.index)
        
    df = df[df.id!='id']
    df = df.drop_duplicates(keep='last')
    df = df.reset_index(drop=True)

    return df

#df = getLevel2info()
#print df.head(20)
