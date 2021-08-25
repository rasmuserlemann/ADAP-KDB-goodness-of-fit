import pandas as pd
import json
import csv

d = pd.read_csv("/Users/rerleman/Dropbox/My Mac (CCI00BHV2JALT)/Documents/distributions_disease.txt", sep='\t', nrows=50)

def sampletransform(x):
    csample = []
    dbsample = []
    for key in x:
        csample = csample + [key]*x[key]['clusterValue']
        dbsample = dbsample + [key]*x[key]['dbValue']
    return([csample,dbsample])

PreviousId = 0
Id = 0
d2 = []
for index, row in d.iterrows():
    samples = sampletransform(json.loads(row['Distribution']))
    if len(samples[0])>2:
        d2.append((samples[0]))
    d.loc[index,'Distribution'] = ','.join(samples[0])

with open("out.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(d2)

#d2.to_csv("/Users/rerleman/Dropbox/My Mac (CCI00BHV2JALT)/Documents/distributions.v4withpvalue.tsv", sep = '\t')

