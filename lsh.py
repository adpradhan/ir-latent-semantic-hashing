import pandas as pd
import numpy as np
import random
import time
import pickle
import copy
# SHINGLE MATRIX CONSTRUCTION
def getShingleMatrix(data,shingle_size):
  SHINGLE_SIZE = shingle_size
  shingle_dict = {}
  cnt = 0
  for seq in data:
      for i in range(0,len(seq)-SHINGLE_SIZE):
          shi = seq[i:i+SHINGLE_SIZE]
          if shi not in shingle_dict:
              shingle_dict[shi] = set()
              cnt += 1
  for j in range(0,len(data)):
    seq = data[j]
    for i in range(0,len(seq)-SHINGLE_SIZE):
        shi = seq[i:i+SHINGLE_SIZE]
        shingle_dict[shi].add(j)
  return shingle_dict

# SIGNATURE MATRIX CONSTRUCTION
def getSignature(permutation, arr, N):
  p_size = len(permutation)
  
  sign = np.full(N,fill_value=p_size+2,dtype=int)
  j = 0
  for shi in arr:
    for i in arr[shi]:
      if permutation[j]<sign[i]:
        sign[i] = permutation[j]
    j += 1 
  return sign

def getSignatureMatrix(shingle_dict,no_of_hashes,N):
  M = len(shingle_dict)
  no_of_hashes = 100
  signature_matrix = []
  number_arr = np.arange(start=1, stop=M+1, step=1)
 
  for i in range(0,no_of_hashes):
    np.random.seed(i)
    hashfn = np.random.permutation(number_arr)
    signature_matrix.append(getSignature(hashfn,shingle_dict,N))
  return signature_matrix

# storing signature_matrix
def storeSignatures(signature_matrix):
  fw = open('signature_policy' , 'wb')
  pickle.dump(signature_matrix, fw)
  fw.close()

# BANDS AND CANDIDATE PAIRS (BUCKETS)
def getBands(b,r,signature_matrix):
  bands = []
  for i in range(b):
    band = signature_matrix[(i*r):((i+1)*r)]
    bands.append(band)
  bands = np.array(bands)
  return bands

def getAdjMatrix(b,r,signature_matrix):
  bands = getBands(b,r,signature_matrix)
  final_bucket = []
  for band in bands:
    k_bucket = {}
    for i in range(len(band[0])):
      bucket_id = hash(tuple(band[:,i]))
      if bucket_id in k_bucket:
        k_bucket[bucket_id].append(i)
      else:
        k_bucket[bucket_id] = [i]
    final_bucket.append(k_bucket)
  
  return final_bucket

# JACCARD SIMILARITY
def jaccard(doc1,doc2,shingle_dict):
  num = 0
  den = 0
  for shi in shingle_dict.keys():
    if doc1 in shingle_dict[shi] or doc2 in shingle_dict[shi]:
      den += 1
    if doc1 in shingle_dict[shi] and doc2 in shingle_dict[shi]:
      num += 1
  if den == 0:
    return 0
  return num/den

# RETRIEVE SIMILAR DNAs TO QUERY
def getSimilarDNA(query,adj,t,shingle_dict):
  ans = []
  for a in adj[query]:
    sim = jaccard(a,query,shingle_dict)
    if sim >= t and a!=query:
      ans.append((a,sim))
  return ans

# For building or loading minhash (signature matrix)
def util(data,t,hashes,shingle_size,train):
  if t == 0.8:
    b = 10
    r = 10
  if t == 0.9:
    b = 5
    r = 20
  if t == 0.55:
    b = 20
    r = 5
  N = len(data)
  
  if train:
    print("Started Shingle Matrix Construction...")
    shingle_dict = getShingleMatrix(data,shingle_size)
    fw = open('shingle_policy' , 'wb')
    pickle.dump(shingle_dict, fw)
    fw.close()
  else:
    fr = open('shingle_policy', 'rb')
    shingle_dict = pickle.load(fr)
    fr.close()
  
  if train:
    print('Started Signature Matrix Construction...')
    signature_matrix = getSignatureMatrix(shingle_dict,hashes,N)
    storeSignatures(signature_matrix)
  else:
    # load signature_matrix
    fr = open('signature_policy', 'rb')
    signature_matrix = pickle.load(fr)
    fr.close()

  if train:
    print('Started Buckets Construction...')
    final_bucket = getAdjMatrix(b,r,signature_matrix)
    fw = open('adj_policy' , 'wb')
    pickle.dump(final_bucket, fw)
    fw.close()
  else:
    fr = open('adj_policy', 'rb')
    final_bucket = pickle.load(fr)
    fr.close()
  
  return (shingle_dict,signature_matrix,final_bucket)

# Fetches similar documents for a given query
def getQueryResult(data,query,t,hashes,shingle_size,train):
  if t == 0.8:
    b = 10
    r = 10
  if t == 0.9:
    b = 5
    r = 20
  if t == 0.55:
    b = 20
    r = 5
  
  query_no = len(data)
  
  (shingle_dict,signature_matrix,final_bucket) = util(data,t,hashes,shingle_size,train)
  
  for i in range(0,len(query)-shingle_size):
    s = query[i:i+shingle_size]
    if s in shingle_dict:
      shingle_dict[s].add(len(data))

  M = len(shingle_dict)
  no_of_hashes = 100
  number_arr = np.arange(start=1, stop=M+1, step=1)
  query_sign = []
  for i in range(0,no_of_hashes):
    np.random.seed(i)
    hashfn = np.random.permutation(number_arr)
    val = len(number_arr) + 2
    j = 0
    for i in shingle_dict:
      if query_no in shingle_dict[i]:
        if val > hashfn[j]:
          val = hashfn[j]
      j+=1
    # signature_matrix[i].append(val)
    query_sign.append(val)

  # print(len(final_bucket))
  bno = 0
  results = set()
  for i in range(0,b):
    # bno = int(i/r)
    hashVal = hash(tuple(query_sign[i*r:(i+1)*r]))
    if hashVal in final_bucket[i]:
      for z in final_bucket[i][hashVal]:
        results.add(z)
     
  final_result = []
  for i in results:
    x = jaccard(query_no,i,shingle_dict)
    if x>=t:
      final_result.append((i,x))

  return final_result


    

	