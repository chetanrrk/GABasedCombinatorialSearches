#!/usr/bin/python
'''
Created on May 29, 2011

@author: Chetan Raj Rupakheti
'''
import numpy as np
import scipy.spatial.distance as sd
import scipy as sc
import time

class MaxMin():
    '''
    Constructs MaxiMin object with given parameters.
    
    @param pool: A 2-D array containing integer sequences that represents a pool of compounds.
    @param size: The selection size from the pool to be computed using MaxiMin algorithm.
    '''
    def __init__(self, pool, size,distType):
        self._pool = pool
        self._size = size
        self._distType = distType
        
    def getMinVal(self,j,distMat,selection):
        minVal=10000000000000
        for key in selection.keys():
            if distMat[key][j]<=minVal:
                minVal = distMat[key][j] 
        return minVal
  
    def getMaxminSet(self):
        selection={}
        print "Computing Distance"
        distMat = sd.cdist(self._pool,self._pool,self._distType)
        print "Starting Maxmin"
        selection[0]=self._pool[0] # append first compound from pool to the selection
        maxVal,maxInd = np.max(distMat[0]),np.argmax(distMat[0])
        selection[maxInd] = self._pool[maxInd]
        minimums = {}
        for i in range(2,self._size):
            for j in range(len(self._pool)):
                minVal = self.getMinVal(j, distMat, selection) # get minimum values for compounds in selection 
                minimums[j]= minVal
            maxVal = np.max(minimums.values()) # compute max of the minimums
            maxIndex = np.argmax(minimums.values()) # get index of maximum compounds
            maxIndex = minimums.keys()[maxIndex] # get actual index relating to the pool
            maxCompound = self._pool[maxIndex] # get maximum compound from pool using proper index
            selection[maxIndex] = maxCompound # append in selection both index and compound
            minimums = {}
        print "Maxmin Done"    
        return selection.values()               
    
        
# Driver for the MaxiMin algorithm with sample test cases

if __name__ == "__main__":
    startTime=time.time()
#    array = [[0,0,0,0],[1,1,0,0],[1,0,1,0],[1,1,1,1],[1,1,1,0],[0,0,0,1]]
    #array = [[1.0, 1.0, 0.0, 0.0, 0.0, 1.0], [1.0, 0.0, 1.0, 1.0, 1.0, 0.0], [1.0, 1.0, 1.0, 0.0, 0.0, 1.0], [1.0, 0.0, 1.0, 0.0, 1.0, 0.0], [1.0, 1.0, 0.0, 1.0, 1.0, 0.0], [0.0, 0.0, 1.0, 1.0, 1.0, 0.0], [1.0, 1.0, 0.0, 0.0, 1.0, 1.0], [0.0, 1.0, 0.0, 0.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [1.0, 0.0, 0.0, 1.0, 0.0, 1.0], [1.0, 0.0, 1.0, 0.0, 1.0, 1.0], [1.0, 0.0, 1.0, 0.0, 1.0, 0.0], [0.0, 1.0, 1.0, 0.0, 1.0, 0.0], [1.0, 1.0, 0.0, 1.0, 0.0, 1.0], [1.0, 0.0, 0.0, 1.0, 0.0, 1.0], [0.0, 1.0, 1.0, 1.0, 0.0, 1.0], [0.0, 0.0, 1.0, 0.0, 1.0, 1.0], [1.0, 0.0, 1.0, 0.0, 1.0, 0.0], [1.0, 1.0, 1.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0, 1.0, 1.0], [1.0, 1.0, 0.0, 1.0, 0.0, 0.0], [1.0, 0.0, 0.0, 1.0, 1.0, 1.0], [0.0, 1.0, 0.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0, 0.0, 0.0], [0.0, 1.0, 1.0, 1.0, 1.0, 1.0], [0.0, 1.0, 0.0, 0.0, 1.0, 1.0], [0.0, 1.0, 1.0, 0.0, 0.0, 1.0], [1.0, 1.0, 1.0, 0.0, 1.0, 0.0], [1.0, 1.0, 0.0, 1.0, 0.0, 1.0], [1.0, 0.0, 1.0, 0.0, 0.0, 1.0], [0.0, 0.0, 1.0, 1.0, 0.0, 1.0], [1.0, 0.0, 0.0, 1.0, 0.0, 1.0], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 0.0, 0.0, 1.0, 1.0], [1.0, 0.0, 0.0, 1.0, 0.0, 1.0], [1.0, 0.0, 0.0, 1.0, 1.0, 0.0]]

#    array = [1,2,3,4,5,6,7,8,9,10]
    #array = [[1,1,1,1],[0,0,0,0],[1,1,0,0],[1,0,1,0],[1,1,1,0],[0,0,0,1]]
#    array = [[0.,  0.,  1.,  0.,  0.,  1.,  1.,  1.,  0.,  0.,  1.,  0.,  0., 0.,  0.,  0.,  0.,  1.,  0.,  0.], [ 1.,  0.,  1.,  1.,  1.,  0.,  0.,  0.,  1.,  0.,  1.,  0.,  1.,
#        0.,  1.,  0.,  0.,  1.,  0.,  0.], [ 1.,  0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  1.,  0.,  1.,  1., 0.,  0.,  0.,  1.,  1.,  0.,  0.], [ 0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  0., 1.,  1.,  0.,  0.,  0.,  0.,  0.], [ 1.,  1.,  0.,  1.,  1.,  1.,  0.,  1.,  1.,  1.,  0.,  1.,  0.,
#        0.,  1.,  0.,  0.,  0.,  0.,  1.], [ 0.,  0.,  1.,  0.,  1.,  1.,  1.,  0.,  1.,  1.,  0.,  1.,  0.,
#        0.,  1.,  0.,  0.,  0.,  0.,  0.]]
    #array = [[0, 0, 0, 0, 0], [0, 0, 0, 0, 1], [0, 0, 0, 1, 0], [0, 0, 0, 1, 1], [0, 0, 1, 0, 0], [0, 0, 1, 0, 1], [0, 0, 1, 1, 0], [0, 0, 1, 1, 1], [0, 1, 0, 0, 0], [0, 1, 0, 0, 1], [0, 1, 0, 1, 0], [0, 1, 0, 1, 1], [0, 1, 1, 0, 0] ,[0, 1, 1, 0, 1], [0, 1, 1, 1, 0], [0, 1, 1, 1, 1], [1, 0, 0, 0, 0], [1, 0, 0, 0, 1] ,[1, 0, 0, 1, 0] ,[1, 0, 0, 1, 1] ,[1, 0, 1, 0, 0], [1, 0, 1, 0, 1], [1, 0, 1, 1, 0], [1, 0, 1, 1, 1] ,[1, 1, 0, 0, 0], [1, 1, 0, 0, 1], [1, 1, 0, 1, 0], [1, 1, 0, 1, 1], [1, 1, 1, 0, 0], [1, 1, 1, 0, 1], [1, 1, 1, 1, 0], [1, 1, 1, 1, 1]]
    array = np.random.randint(0,2,size = (1000,20))
    print np.shape(array)    
    maxmin = MaxMin(array,100,'Hamming')
    #maxmin.compute()
    #print maxmin.compute()
    print len(maxmin.getMaxminSet())
    print "total time",time.time()-startTime    