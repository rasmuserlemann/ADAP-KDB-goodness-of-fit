import pandas as pd
import json
from scipy.stats import ks_2samp

d = pd.read_csv("/Users/rerleman/Dropbox/My Mac (CCI00BHV2JALT)/Documents/distributions.tsv", sep='\t', nrows=50)
d['KSpvalue'] = ''

def sampletransform(x):
    csample = []
    dbsample = []
    for key in x:
        csample = csample + [key]*x[key]['clusterValue']
        dbsample = dbsample + [key]*x[key]['dbValue']
    return([csample,dbsample])

PreviousId = 0
for index, row in d.iterrows():
    samples = sampletransform(json.loads(row['Distribution']))
    d.loc[index,'KSpvalue'] = ks_2samp(samples[0],samples[1])[1]
    
d.to_csv("/Users/rerleman/Dropbox/My Mac (CCI00BHV2JALT)/Documents/distributions.v4withpvalue.tsv", sep = '\t')

