from random import randint
import numpy as np
import yaml 

'''
The coordinates of cell positions will be calculated in 'angstrom' and 
will be converted into 'nanometer' .

This is because of randint (random integer).
'''

def genCellCoord2D(NoCell=10,
                   Rlim=0.1,
                   lowx=0,
                   highx=10,
                   lowy=0,
                   highy=10,
                   outputName='test',
                   restingRatio=0.5):
    
    numP = 0
    maxP = int(NoCell)
    maxIter = 1000*maxP
    minDist = Rlim
    loopcounter = 1
    NoOfResting = int(NoCell*restingRatio)
    coordlib = {}
    xo = np.zeros(maxP)
    yo = np.zeros(maxP)
    
    # recording a number of cells in the simulation box 
    coordlib['NoCell'] = NoCell
    
    while numP < maxP and loopcounter < maxIter:
        xpossible = randint(lowx+1,highx-1)
        ypossible = randint(lowy+1,highy-1)
        if numP == 0:
            xo[numP] = xpossible
            yo[numP] = ypossible
            marker = 'resting'
            cellNo = str(numP)
            coordlib[cellNo] = [xpossible/highx,ypossible/highy,marker]
            numP = numP + 1
            continue

        distance = np.sqrt((xo-xpossible)**2 + (yo-ypossible)**2)
        
        if min(distance) >= minDist:
            xo[numP] = xpossible
            yo[numP] = ypossible
            
            #--- assigning state of microlgia ---#
            if numP <= NoOfResting:
                marker = 'resting'
            elif numP > NoOfResting: 
                marker = 'activated'
            #------------------------------------#
            
            cellNo = str(numP)
            coordlib[cellNo] = [xpossible/highx,ypossible/highy,marker]
            numP = numP + 1
        loopcounter = loopcounter + 1
        
    with open(outputName+'.yml','w') as file:
      document = yaml.dump(coordlib,file)

    return 

#!/usr/bin/env python
import sys
import numpy as np
##################################
#
# Revisions
#       10.08.10 inception
#
##################################
#
# Message printed when program run without arguments
#
def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose: 

"""
  return msg

#
# MAIN routine executed when launching this script from command line
#
if __name__ == "__main__":
  import sys
  msg = helpmsg()
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  NoCell=10
  Rlim=5
  lowx=0
  highx=100
  lowy=0
  highy=100
  outputName='test'
  restingRatio=0.5
    
  for i,arg in enumerate(sys.argv):
    # calls 'runParams' with the next argument following the argument '-validation'
    if arg=="-NoCell":
        NoCell = np.int(sys.argv[i+1])
        
    if arg=="-outputName":
        outputName = sys.argv[i+1]
    
    if arg=="-restingRatio":
        restingRatio = np.float(sys.argv[i+1])
         
genCellCoord2D(NoCell=NoCell,
               Rlim=Rlim,
               lowx=lowx,
               highx=highx,
               lowy=lowy,
               highy=highy,
               outputName=outputName,
               restingRatio=restingRatio
               )