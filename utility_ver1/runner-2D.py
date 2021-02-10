import os
import yaml
import numpy as np
import gc

'''
The order of calculation in single run
[Test mode]
1. Calculate the external ATP signal (fenics_origin-2D)
2. Generate Coordinates of Cell positions 
3. Calculate cellular ATP signals as function of external ATP signals (applied in boundary condition)
'''
########################################
#                                      #
#    Control Panel                     #
#                                      #
########################################
Steps1 = 3                                  # number of steps within 2 ms simulation with stimulation on
Steps2 = 10                                 # number of steps within 498 ms simulation with no stimulation
Duration1 = 0.002                           # simulation length in second for stimulation-on calculation
Duration2 = 0.498                           # simulation length in second for stimulation-off calculation
boxLength = 500                             # box length in um
NoCell = 20                                 # number of cells in the box 
total_jobLength = 1200                         # total job length by steps 
jobcontinued = 0                            # job you desired continued from ___ This is very latest test.yml number
jobLength = total_jobLength - jobcontinued  # actual number of calculation steps 
autocrine = 0                               # autocrine on with other than 0 or off with 0
testmode = 1
########################################

## Extract the simulation condition 
path = '../utility/'
singularity = 'singularity exec /home/bending456/singularity-img/fenics_ben_2019.img python3 '

# step 1: generate the coordinates of cells
if jobcontinued == 0:
    if testmode == 0:
        os.system('python3 '+path+'coordGen-2D.py -outputName test0 -NoCell '+ str(NoCell))
    else:
        os.system('cp ../inputs/test0.yml ./')
    print('---- step 1 completed ----')
    # step 2: calculation of external ATP signal 
    os.system(singularity+path+'fenics-2D-refined.py -InputCoords test0 -T '+str(Duration1)+' -steps '+str(Steps1)+' -CalcOut run1 -lengthOfBox '+str(boxLength)+' -stim 1 -autocrine '+str(autocrine))
    print('---- step 2 completed ----')
    print('---- the initial 2 steps are completed ----')

n = 2
stateVar_initial = {}
for c in np.arange(NoCell):
    stateVar_initial[str(c)] = 0
    
with open('stateVar0.yml','w') as file:
        document = yaml.dump(stateVar_initial,file) 

for i in np.arange(jobLength):  
    # step 3: displacement 
    n = n + 1
    
    if jobcontinued != 0:
        i = i + jobcontinued
        
    if (i+1)%2 == 0:
        # Turn on DBC
        stim = 1
        Duration = Duration1
        DispDelta_t = Duration1*(i+1)+Duration2*i
        Steps = Steps1
        print('------- stimulation on -----------')
    else:
        # Turn on NBC
        stim = 0
        Duration = Duration2
        DispDelta_t = Duration1*(i+1)+Duration2*(i+1)
        Steps = Steps2
        print('------- stimulation off ----------')
    
    os.system(singularity+path+'displacement.py -InputCoords test'+str(i)+' -OutputCoords test'+str(i+1)+' -InputU run'+str(i+1)+' -InputMesh mesh -time '+str(DispDelta_t)+' -intvl '+str(Duration)+' -InputStateVar stateVar'+str(i)+' -OutputStateVar stateVar'+str(i+1)+' -lengthOfBox '+str(boxLength))    
    print('---- step %d: Displacement process is completed ----' %(n))
    gc.collect()
        
    os.system(singularity+path+'fenics-2D-refined.py -InputCoords test'+str(i+1)+' -T '+str(Duration)+' -steps '+str(Steps)+' -CalcOut run'+str(i+2)+' -PrevCalc run'+str(i+1)+' -lengthOfBox '+str(boxLength)+' -stim '+str(stim)+' -Auto '+str(autocrine))
    print('---- step %d: Diffusion calculation is completed ----' %(n))
    gc.collect()
        
#import subprocess
#name = ""
#for i in np.arange(Steps):
#    name += 'run'+str(i+1)+'.png '

#subprocess.call("convert -delay 25 -loop 5 " + name +" " + GifOutPut, shell=True)
#os.system('rm -rf *.png')
#os.system('rm -rf *.h5')
#os.system('rm -rf *.yml')
