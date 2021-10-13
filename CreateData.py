import pandas as pd
import json
import csv

d = pd.read_csv("/Users/rerleman/Documents/adap-kdp-gof-data/data.tsv", sep='\t')

def sampletransform(x):
    csample = []
    dbsample = []
    for key in x:
        csample = csample + [key]*x[key]['clusterValue']
        dbsample = dbsample + [key]*x[key]['dbValue']
    return([csample,dbsample])

PreviousId = 0
Id = 0
d2 = d
d2["Cluster"] = ""
d2["DB"] = ""
d2["P-value"] = ""
for index, row in d.iterrows():
    samples = sampletransform(json.loads(row['Distribution']))
    d2.loc[index,'Cluster'] = ','.join(samples[0])
    d2.loc[index,'DB'] = ','.join(samples[1])
    if index % 10000 == 0:
        print(index)

d2.to_csv("/Users/rerleman/Documents/adap-kdp-gof-data/data1.tsv", sep = '\t')

