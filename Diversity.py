'''
Created on Jan 28, 2014

@author: chetan
'''
import os,sys
import numpy as np
import scipy as sc
import scipy.spatial.distance as sd
from MaxMin import *
from heapq import *
import copy
import itertools
import generateCombination
import pickle

nkLandScape = generateCombination.nkLandScape


'''
NK fitness landscape variables
'''    
k = 9 ## links of a cell to other cells
N = 19 ## length of the binary array

'''
GA Run variables
'''    
fType = 'nk'  ## choice of the fitness landscape
dType = 'hamming' ## distance metric
trials = 1 ## number of GA trials
generations = 50 ## number of GA iterations
xoverP = 0.10 ## probability of cross over
mutP = 0.05   ## probability of mutation
numPos = 4 #number of positions to flip the bit during mutation
minSize = 100 #minimum number of solutions that must be present every generation
SubSetSize = 200 ## number of solutions to retain each iteration


'''
Initial Solutions
'''
lowerLimit=0;upperLimit=2;collectionSize=200;vec_dim=N


'''
Generates fitness table to read fitness values from
@param k: number of links of a cell to other cells
'''
def generateFitnessTable(k):
    comb = makeComb(k)
    fitnessTable={}
    for c in comb:
        fitnessTable[c] = np.random.rand()
    return fitnessTable

'''
computes fitness of a binary array using griewank fitness function
@param vec: binary array
'''
def griewank(vec):
    temp1=0
    temp2=1
    for i in range(len(vec)):
        temp1 = temp1+((vec[i]-100)*(vec[i]-100))
        temp2 = temp2*(np.cos((vec[i]-100)/np.sqrt(i+1))+1)
    Fx = (temp1-temp2)/4000;
    return Fx

'''
computes fitness of a binary array using rastrigin fitness function
@param vec: binary array
'''
def rastrigin(vec):
    Fx = 0
    for i in range(len(vec)):
        temp = (vec[i]*vec[i])-(10*np.cos(2*np.pi*vec[i]))+10
        Fx=Fx+temp
    return Fx 

'''
generates initial guess of the binary array solutions
@param lowerLimit: 0; upperLimit:1
@param collectionSize: number of initial solutions at the start
@param vec_dim: length of the binary array
'''
def makeVec(lowerLimit,upperLimit,collectionSize,vec_dim):
    X = np.random.randint(lowerLimit,upperLimit,size = (collectionSize,vec_dim))
    return X    

'''
Selects fittest solution by roulette wheel selection
@param vecSet: collection of all solutions
@param fType: fitness function type (eg nk)
@param k: number of link of a cell to other cells
@param fitnessTable: table that contains fitness value of an array
'''
def rouletteWheel(vecSet,fType,k,fitnessTable):
    fitnessArr = []
    selected = {}
    i = 0 # index for selected set
    j = 0 # index for vector set
    for vec in vecSet:
        if fType == 'gr':
            fitnessArr.append(griewank(vec))
        elif fType == 'rs':
            fitnessArr.append(rastrigin(vec))
        elif fType == 'nk':
            fitnessArr.append(nkLandScape(vec,k,fitnessTable))    
    for fitness in fitnessArr:
        scaledFitness = (fitness-np.mean(fitnessArr))/np.std(fitnessArr)
        if scaledFitness >= np.random.rand():
            selected[i]=vecSet[j],fitness
            i+=1
        j+=1
    return selected

def tournamentSelection(vecSet,selectionSize,fType,k,fitnessTable):
    selected = {}
    for i in range(selectionSize):
        num1 = np.random.randint(0,len(vecSet))
        num2 = np.random.randint(0,len(vecSet))
        if fType =='gr':
            if griewank(vecSet[num1]) > griewank(vecSet[num2]):
                selected[i] = vecSet[num1], griewank(vecSet[num1])
            else:
                selected[i] = vecSet[num2], griewank(vecSet[num2])
        elif fType =='rs':        
            if rastrigin(vecSet[num1]) > rastrigin(vecSet[num2]):
                selected[i] = vecSet[num1], rastrigin(vecSet[num1])
            else:
                selected[i] = vecSet[num2], rastrigin(vecSet[num2])
        elif fType =='nk':
            if nkLandScape(vecSet[num1],k,fitnessTable) > nkLandScape(vecSet[num2],k,fitnessTable):
                selected[i] = vecSet[num1], nkLandScape(vecSet[num1],k,fitnessTable)
            else:
                selected[i] = vecSet[num2], nkLandScape(vecSet[num2],k,fitnessTable)             
    return selected                

def distance(arr1, arr2,dtype):
    '''
    Returns the distance between guess and code array.
    '''
    if len(arr1)>1 and len(arr2)>1:
        return sd.cdist(arr1, arr2,dtype)
    else:        
        return sd.cdist([arr1],[arr2],dtype)[0][0]
        
    
def nearestNeighborDist(mySet,dType):
    dMatX=sd.cdist(mySet,mySet,dType)
    minD=[]
    j=0
    for i in range(len(dMatX)):
        arr=dMatX[i]
        ind = sc.argmin(arr)
        if ind == j:
            arr = np.delete(arr,ind)
            myMin,ind = np.min(arr),sc.argmin(arr)
            minD.append(myMin) 
        j+=1
    nnDist = float(np.sum(minD)/len(minD))
    return nnDist        

    #older and slower implementation of nearest-neighbor distance          
    """            
    nnDist=0
    min=[] # stores nearest neighbor distance for each compound 
    minD = 100000 
    for i in range(len(mySet)):
        for j in range (len(mySet)):
            if i!=j:
                if dType=='Ham':
                    d = distance(mySet[i],mySet[j])
                elif dType == 'Euc':
                    d = eucDist(mySet[i],mySet[j])    
                '''
                Compute running minimum
                '''
                if d<minD:
                    minD = d
        min.append(minD)
        minD =  100000       
    nnDist = np.sum(min)
    #print len(mySet)
    return float(nnDist/len(mySet))
    """    

def getVec(vec):
    myVec=[]
    for ele in vec:
        myVec.append(ele[0])
    return myVec     


def selectFit(threshold,vecSet,fType,k,fitnessTable):
    wholeSet=[]
    selectedArr=[]
    for ele in vecSet:
        if fType == 'gr':
            wholeSet.append((griewank(ele),ele))
        elif fType == 'rs':
            wholeSet.append((rastrigin(ele),ele))
        elif fType == 'nk':
            wholeSet.append((nkLandScape(ele,k,fitnessTable),ele))    

    for ele in wholeSet:
        if ele[0]>=threshold:
            selectedArr.append(ele[1])    
    return selectedArr                
            
def averageFitness(vecSet,fType,k,fitnessTable):
    fitnessArr=[]
    for vec in vecSet:
        if fType=='rs':
            fitnessArr.append(rastrigin(vec))
        elif fType=='gr':
            fitnessArr.append(griewank(vec))
        elif fType=='nk':
            fitnessArr.append(nkLandScape(vec,k,fitnessTable))            
    fitness = np.sum(fitnessArr)
    return fitness/len(vecSet)            

def maxFit(vecSet,fType,k,fitnessTable):
    fitnessArr=[]
    for vec in vecSet:
        if fType=='rs':
            fitnessArr.append(rastrigin(vec))
        elif fType=='gr':
            fitnessArr.append(griewank(vec))
        elif fType=='nk':
            fitnessArr.append(nkLandScape(vec,k,fitnessTable))            

    return np.max(fitnessArr),vecSet[np.argmax(fitnessArr)]          



'''
This function performs single point crossover where the position to do cross over is determined randomly
'''
def singlePointCrossOver(parentSet):    
    #print parentSet
    vec1 = parentSet[0]
    vec2 = parentSet[1]
    a = np.random.randint(1,len(vec1)-1) # determines where the split will occur
    c = np.concatenate((vec1[0:a],vec2[a:len(vec2)]),axis=0)#vec1[0:a]+vec2[a:len(vec2)]
    d = np.concatenate((vec2[0:a],vec1[a:len(vec1)]),axis=0)       
    Children = [c,d]
        
    return Children


def doMutation(mutP):    
    if mutP > np.random.rand():
        return True
    else:
        return False
'''
def Mutate(set,mutP,numPos):

    mutants = []
    for guess in set:     
        """
        we have to try to mutate 'numPos' bit in a chromosome as allowed by mutation probability
        """
        for i in range (numPos):
            pos = np.random.randint(0,len(guess)-1)
            if doMutation(mutP):    
                if guess[pos] == 0:
                    guess[pos] = 1
                elif guess[pos] == 1:
                    guess[pos] = 0
        mutants.append(guess)    
                 
    return mutants
'''

def Mutate(parents,mutP,numPos):
    '''
    return the parent and the new children in the end not only the mutants? 
    '''
    mutants = []
    for guess in parents:     
        '''
        we have to try to mutate 'numPos' bit in a chromosome as allowed by mutation probability
        '''
        for i in range (numPos):
            pos = np.random.randint(0,len(guess)-1)
            if doMutation(mutP):    
                if guess[pos] == 0:
                    guess[pos] = 1
                elif guess[pos] == 1:
                    guess[pos] = 0
        mutants.append(guess)    
    children = mutants+parents         
    return children

def doCrossOver(xOverP):
    if xOverP > np.random.rand():
        return True
    else:
        return False
    
def getSolution(selectedChildren):
    solutions=[]
    for c in selectedChildren.values():
        solutions.append(c[0])   
    return solutions    

'''
How do you set the threshold???How do you compare??
'''
        
def setThreshold(thresholdOld,mySet,k,fitnessTable,generations,minSize):
    target= 0.736483402261 #lets make this the max value
    fitnessArr=[]
    for ele in mySet:
        fitnessArr.append(nkLandScape(ele,k,fitnessTable))
    mean = np.mean(fitnessArr)
    stdev = 2*np.std(fitnessArr) # lets weedout the tails 4.4% assuming gaussian
    
    threshold = thresholdOld
    threshold = mean+mean*0.01 
    if threshold >= target:
        threshold = target
    '''    
    else:
        print 'Increasing Threshold...old threshold was',threshold
        threshold = mean+mean*0.01  
        print 'New Threshold is',threshold  
    #elif generations>len(mySet[0])/2 and len(mySet)>minSize:
    #    threshold = thresholdOld+thresholdOld*0.01    
    '''    
    return threshold

def copyParents(parents,vecSet,minSize,fType):
    if len(vecSet)<minSize:
        for i in range(minSize-len(vecSet)):
                
            val = np.random.randint(0,len(parents))
            #print "index of parent picked is",val
            vecSet.append(parents[val])
    return vecSet

def makeComb(num):
    return ["".join(seq) for seq in itertools.product("01", repeat=num)]
        
def stringArr(arr):
    tempi = []
    tempf = []
    for ele1 in arr:
        for e in ele1:
            tempi.append(e)    
        tempf.append(tempi)
        tempi=[]
    return tempf

def readFitnessTable(myFile,N,k):
    #fin=open('N'+str(N)+'K'+str(k)+'_FitnessTable.txt','r')
    fin = open(myFile,'r')
    fitnessTable={}
    for line in fin:
        token=line.rstrip().split()
        key=token[0];val=float(token[1])
        fitnessTable[key]=val
    return fitnessTable    

def increaseMutationRate(mutP,N,j,vecSet,minSize):
    
    if j<N/2 or len(vecSet)<minSize:
        newMutP = mutP #+2*mutP
        return newMutP
    else:
        print "Mutation Here"
        newMutP = 0.15
        return newMutP

def pickXOver(vecset,num):
    pairs=[]
    for i in range(num):
        pairs.append(vecSet[np.random.randint(len(vecSet))])
    return pairs
        
def run(pool,vecSet,mutP,selectionLength):
    
    children=[]
    for i in range(len(vecSet)):
        pairs = pickXOver(vecSet,2) #selecting two parents for x-over  
        if doCrossOver(xoverP):
            tmp = singlePointCrossOver(pairs)
            children.append(tmp[0]);children.append(tmp[1])
        else:    
            tmp = copy.deepcopy(pairs)
            children.append(tmp[0]);children.append(tmp[1])
                
    mutants = Mutate(children,mutP,numPos)
    
    for c in mutants:
        try:
            c = c.tolist()
        except AttributeError:continue              
        if c not in pool:
            pool.append(c)
    #print "length of uniqC is",len(uniqC)        

    if len(pool)>selectionLength:            
        maxmin = MaxMin(pool,selectionLength,dType)
        selectedChildren = maxmin.getMaxminSet()
        return selectedChildren,pool
    else:
        return pool,pool    

def stdFitness(vecSet,fType,k,fitnessTable):
    fitnessArr=[]
    for vec in vecSet:
        if fType=='rs':
            fitnessArr.append(rastrigin(vec))
        elif fType=='gr':
            fitnessArr.append(griewank(vec))
        elif fType=='nk':
            fitnessArr.append(nkLandScape(vec,k,fitnessTable))            
    stdFitness = np.std(fitnessArr)
    return stdFitness


if __name__=="__main__":
    restart=True
    pool = []
    fout = open('Diversity1.txt','w')
    '''
    Initial Solutions
    '''
    lowerLimit=0;upperLimit=2;collectionSize=200;vec_dim=N
    
    vecSet = makeVec(lowerLimit,upperLimit,collectionSize,vec_dim)
    pool = copy.deepcopy(vecSet.tolist())
    selectionLength = 200
    print 'Vector set done!!!';
    print>>fout,'Generation','Diversity','PoolLength','Subset'
     
    if restart == True:
        try:
            pool=pickle.load(open("pool.p","rb"))
            print "Continuing old run"
            lastIter=0
            for file in [doc for doc in os.listdir('.') if doc.startswith("vecSet")]:
                if int(file.split('vecSet')[1].split('.')[0]) > lastIter:
                    lastIter = int(file.split('vecSet')[1].split('.')[0])
                 
            vecSet = pickle.load(open("vecSet"+str(lastIter)+'.p','rb'))
        except IOError:print
            
    for i in range(trials):    
        for j in range(generations):
            print "-------Running generation ",str(j),"-------"
            
            print>>fout,j,nearestNeighborDist(vecSet,dType),len(pool),len(vecSet)
            
            vecSet,pool = run(pool,vecSet,mutP,selectionLength)
            #print 'pool',pool;print 'vecset',vecSet                         
            #pool = pool+vecSet
            print "length of pool",len(pool) 
            if j%5==0:
                pickle.dump(pool,open('pool.p','wb'))
                pickle.dump(vecSet,open('vecSet'+str(j)+'.p','wb'))
            #mutP = increaseMutationRate(mutP,N,j,vecSet,minSize)  
        print "----------------Final Result----------------"
        print>>fout, j,nearestNeighborDist(vecSet,dType),len(pool),len(vecSet)
                
        print "Done!"
        pickle.dump(vecSet,open('diverseSet.p','wb'))
        pickle.dump(pool,open('pool.p','wb'))
        fout.close()

    def printFitness(vecs,fType,fileName,k,fitnessTable):
        fout = open(fileName,'w')
        for ele in vecs:
            if fType=='rs':
                print>>fout, rastrigin(ele)
            elif fType=='gr':    
                print>>fout, griewank(ele)
            elif fType=='nk':
                print>>fout,nkLandScape(ele,k,fitnessTable)    
        fout.close()
