from fenics import *
from mshr import *
import numpy as np 
import yaml
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import curve_fit
from scipy.integrate import odeint
import math

'''
The unit of length in this calculation is in [nanometer]!!
'''

def stateFunc(y,t,c):
    ## Rate of state transition 
    kf = 1
    kb = 0.05
    
    dydt = kf*c*(1-y) - kb*y
    
    return dydt

def displacement(InputCoords='test',
                 OutputCoords='test',
                 InputU='test',
                 InputMesh='test',
                 lengthOfBox='test',
                 time=2,
                 interval=2,
                 InputStateVar='stateVar',
                 OutputStateVar='stateVar'
                ):
    
    
    ## Some constant 
    SearchingR = 5
    
    ## Open the previous coordinates of cells 
    with open(InputCoords+'.yml') as file1:
        coords = yaml.load(file1,Loader=yaml.FullLoader)
    
    NoOfCells = coords['NoCell']
    ## Loading up the previously calculated concentrations and projected onto the previous mesh
    mesh = Mesh(InputMesh+'.xml')        
    uFile = HDF5File(MPI.comm_world,InputU+'.h5','r')
    V = FunctionSpace(mesh,'P',1)
    u = Function(V)
    uFile.read(u,'/u')
    uFile.close()
    u.set_allow_extrapolation(True)
    
    ## Setting up the odeint calculation
    timeBefore = 0 #time - interval
    timeTarget = time 
    t = scipy.linspace(timeBefore,timeTarget,10)
    ## Setting state variable input
    with open(InputStateVar+'.yml') as file2:
        stateVariable = yaml.load(file2,Loader=yaml.FullLoader)

    newStateVar = {}
    
    ## Loading up the local concentration 
    with open('LocalConcOfCell.yml') as file3:
        LocalConc = yaml.load(file3,Loader=yaml.FullLoader)
    
    print('-- State Variable Dependent Factor --')
    print('-- Box Length: '+str(lengthOfBox)+' um --')
    for i in np.arange(NoOfCells):
        xi = coords[str(i)][0]*lengthOfBox
        yi = coords[str(i)][1]*lengthOfBox
        marker = coords[str(i)][2]
        #ui = u(xi,yi)
        ui = LocalConc[str(i)]
        ## Force calculation 
        xcoord = [xi-SearchingR, xi+SearchingR]
        ycoord = [yi-SearchingR, yi+SearchingR]
        # Extracting the concentration information from the previous calculation
        # within the searching parameters
        ##             u12
        ##      u01    u11    u21
        ##             u10
               
        #stateVarc = odeint(stateFunc,stateVariable[str(i)],t,args=(ui,))
        stateVarc = odeint(stateFunc,0,t,args=(ui*(1e18),))
        dummy1 = np.float(stateVarc[-1])
        newStateVar[str(i)] = dummy1
        stateDepFactor = 1/(1+(dummy1/0.65)**5) #1/(1+(dummy1/0.3)**2)
                
        u11 = ui
        
        # Prevent the error at the boundary of simulation box
        if xcoord[0] < 1e-10:
            u01 = 1e-40
        else:
            u01 = u(xcoord[0],yi)
            
        if xcoord[1] > 500 - 1e-14:
            u21 = 1e-40
        else:
            u21 = u(xcoord[1],yi)
            
        if ycoord[0] < 1e-10:
            u12 = 1e-40
        else:
            u12 = u(xi,ycoord[1])
            
        if ycoord[1] > 500 - 1e-14:
            u10 = 1e-40
        else:
            u10 = u(xi,ycoord[0])
        
        ## in x direction 
        ux = (u21 - u01)/2
        ## in y direction 
        uy = (u12 - u10)/2
        
        if abs(ux) > abs(uy):
            dx = ux/abs(ux)
            dy = uy/abs(ux)
        elif abs(ux) <= abs(uy):
            dx = ux/abs(uy)
            dy = uy/abs(uy)
            
        ## displacement
        if marker == 'resting':
            migrationRate = 2e19
        elif marker == 'activated':
            migrationRate = 2e18
            
        newx = np.float(xi+dx*migrationRate*stateDepFactor*ui*interval)
        newy = np.float(yi+dy*migrationRate*stateDepFactor*ui*interval)
        print('conc: '+str(u11)+' stateVar: '+str(dummy1)+' Dependent Factor: '+str(stateDepFactor))
        
        coords[str(i)] = [newx/lengthOfBox,newy/lengthOfBox,marker]
    
    with open(OutputCoords+'.yml','w') as file3:
        document1 = yaml.dump(coords,file3)
       
    with open(OutputStateVar+'.yml','w') as file4:
        document2 = yaml.dump(newStateVar,file4)
        
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

  InputCoords='test'
  OutputCoords='test'
  InputU='test'
  InputMesh='test'
  lengthOfBox=10
  time=2
  interval=2
  InputStateVar='stateVar'
  OutputStateVar='stateVar'
    
  for i,arg in enumerate(sys.argv):
    # calls 'runParams' with the next argument following the argument '-validation'
    if arg=="-InputCoords":
        InputCoords = sys.argv[i+1]
    
    if arg=="-OutputCoords":
        OutputCoords = sys.argv[i+1]
    
    if arg=="-InputU":
        InputU = sys.argv[i+1]
         
    if arg=="-InputMesh":
        InputMesh = sys.argv[i+1]
        
    if arg=="-lengthOfBox":
        lengthOfBox = np.int(sys.argv[i+1])
        
    if arg=="-time":
        time = np.float(sys.argv[i+1])
    
    if arg=="-intvl":
        interval = np.float(sys.argv[i+1])
        
    if arg=="-InputStateVar":
        InputStateVar = sys.argv[i+1]
        
    if arg=="-OutputStateVar":
        OutputStateVar = sys.argv[i+1]
        
displacement(InputCoords=InputCoords,
             OutputCoords=OutputCoords,
             InputU=InputU,
             InputMesh=InputMesh,
             lengthOfBox=lengthOfBox,
             time=time,
             interval=interval,
             InputStateVar=InputStateVar,
             OutputStateVar=OutputStateVar
             )
