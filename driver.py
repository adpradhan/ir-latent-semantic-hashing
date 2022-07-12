import pandas as pd
import time
import numpy as np
from lsh import getShingleMatrix, getSignature, getSignatureMatrix, storeSignatures, getBands, getAdjMatrix, jaccard
from lsh import getSimilarDNA, util, getQueryResult

d = pd.read_csv('human_data.txt',sep="\t")
d = d['sequence']
q = d[119]

t = 0.8
shingle_size = 6
hashes = 100
s = time.time()
x = getQueryResult(d,q,t,hashes,shingle_size,False)
e = time.time()
print(f"Time: {e-s}")
for i in x:
	print(i[0],i[1])