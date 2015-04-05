'''
Created on Jul 25, 2013

@author: chetan
'''
#!/usr/bin/python
import os,sys,gzip
import itertools
import numpy as np
#import Selection
import scipy as sc
import scipy.spatial.distance as sd
import  matplotlib.pyplot as plt


def nearestNeighborDist(mySet,dType,numSlice):
    init=0;final=0
    z=0;minDist=[]
    while final<=1:
        A = mySet
        final = 1.0/numSlice+init
        B = mySet[int(init*len(mySet)):int(final*len(mySet))]
        
        print 'distance matrix started' 
        distMat = sd.cdist(B,A,dType)*len(mySet[0])
        j=0
        for i in range(len(distMat)):
            arr=distMat[i]
            ind = sc.argmin(arr)
            if ind == j:
                arr = np.delete(arr,ind)
                myMin,ind = np.min(arr),sc.argmin(arr)    
            minDist.append(myMin) 
            j+=1 
        init=final
        print "Distance Matrix Computed of Slice ",z
        z+=1 
    nnDist = float(np.sum(minDist)/len(minDist))
    return nnDist 

def generateFitnessTable(k,p,alleleVal):
    comb = makeComb(k)
    fitnessTable={}
    #wts=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
    i=0
    for c in comb:
        fitnessTable[c] = np.random.rand()#wts[i]
        i+=1
    '''
    Lets pick 'p' proportion of fitnessTable to be set to 0
    '''
    unfitLen = int(p*len(fitnessTable))
    print "unfitLen",unfitLen
    i=0    
    for key in fitnessTable:
        if i<unfitLen:
            fitnessTable[key] = alleleVal  
        i+=1     
    return fitnessTable

def GenerateEpiPos(N,k):
    epiStats={}
    temp=[]
    for i in range(N):
        m = i
        #epis= np.random.permutation([a for a in range(N)])
        #temp = epis[:k+1]
        for j in range(k+1):
            l = (m+1)%N
            temp.append(l)
            m=m+1
        
        if i in temp:
            epiStats[i]=temp
        else:
            #temp = temp.tolist()
            temp.remove(temp[len(temp)-1])
            temp.append(i)
            epiStats[i]=temp
        temp=[]        
    return epiStats        

'''
Function based on Barnett, et al, 'Ruggedness and Neutrality- The NKp family of Fitness Landscapes'
'''
def nkLandScape(vec,k,fitnessTable):
    epiPos = GenerateEpiPos(len(vec),k)
    fi=0.0
    temp=[]
    fitness=[]
    for i in range(len(vec)):
        for j in epiPos[i]:
            temp.append(vec[j])
        #print "".join(str(i) for i in temp),fitnessTable["".join(str(i) for i in temp)]       
        fi = fi+fitnessTable["".join(str(i) for i in temp)]
        '''    
        if fType == 'gr':     
            fi = fi+griewank(temp)*np.random.rand()#(1.0/(len(vec)))
            fitness.append(fi)
        elif fType == 'rs':
            fi = fi+rastrigin(temp)*np.random.rand()#(1.0/(len(vec)))    
            fitness.append(fi)
        elif fType == 'cum':
            fi = fi+cumFit(temp)*np.random.rand()#(1.0/(len(vec)))    
        '''
        temp=[]
    return fi/len(vec)#np.sum(fitness)
    
def cumFit(vec):
    fi=0.0
    for v in vec:
        fi=fi+v   
    return fi
    
def griewank(vec):
    temp1=0
    temp2=1
    
    for i in range(len(vec)):
        temp1 = temp1+((vec[i]-100)*(vec[i]-100))
        temp2 = temp2*(np.cos((vec[i]-100)/np.sqrt(i+1))+1)
    Fx = (temp1-temp2)/4000;
    return Fx

def rastrigin(vec):
    Fx = 0
    for i in range(len(vec)):
        temp = (vec[i]*vec[i])-(10*np.cos(2*np.pi*vec[i]))+10
        Fx=Fx+temp
    return Fx 

def makeComb(num):
    return ["".join(seq) for seq in itertools.product("01", repeat=num)]
        
def stringArr(arr):
    tempi = []
    tempf = []
    for ele1 in arr:
        for e in ele1:
            tempi.append(int(e))    
        tempf.append(tempi)
        tempi=[]
    return tempf

def readFitnessTable(N,k):
    fin=open('N'+str(N)+'K'+str(k)+'_FitnessTable.txt','r')
    fitnessTable={}
    for line in fin:
        token=line.rstrip().split()
        key=token[0];val=float(token[1])
        fitnessTable[key]=val
    return fitnessTable    


if __name__=="__main__":
    N=19
    k = 9
    p = 0
    alleleVal = 0 # the fitness value of allele    
                
    print 'Generating Fitness Table'
    #fitnessTable={'00':0.1,'01':0,'10':0,'11':0.5}
    
    fitnessTable = generateFitnessTable(k+1,p,alleleVal)
    print fitnessTable
    fout1=open('N'+str(N)+'K'+str(k)+'P'+str(p)+'2_FitnessTable.txt','w')
    for key in fitnessTable:
        print>>fout1,key,fitnessTable[key]
    fout1.close()      
    
    
    
    fout2=open('N'+str(N)+'K'+str(k)+'P'+str(p)+'2FitnessValues.txt','w')
    print 'Generating Arrays'
    arr1 = makeComb(N)
    arr2 = stringArr(arr1)
    print 'Arrays Generated...Evaluating Fitness'
    for i in range(len(arr2)):
        print>>fout2,arr1[i],nkLandScape(arr2[i],k,fitnessTable)
    sys.exit()
    distMat = sd.cdist(arr2,arr2,'hamming')*N 
    plt.imshow(distMat);plt.colorbar()
    plt.show()

    #sys.exit()
    
    '''
    fitnessTable= readFitnessTable(N,k)
    arr = makeComb(N)
    binArr = stringArr(arr)
    targetDic={}
    i=0
    for ele in arr:
        fitness = nkLandScape(ele,k,fitnessTable)
        if fitness>=0.534:
            targetDic[i]=fitness
            i+=1  
            #print>>fout3,nkLandScape(ele,k,fitnessTable)                  
    targetArr=[]        
    for key in targetDic:
        targetArr.append(binArr[key])
    
    print nearestNeighborDist(targetArr,'Hamming',10)    
    '''       