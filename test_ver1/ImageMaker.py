import subprocess
import numpy as np
import os

nameloc = ""
oloc = 'test-location-auto.gif'

for i in np.arange(1200):
    if (i + 1)%1 == 0:
        nameloc += 'run'+str(i+1)+'.png '

os.system("convert -delay 5 -loop 5 " + nameloc +" " + oloc)

nameall = ""
oall = 'test-all-auto.gif'
for i in np.arange(1200):
    if (i+1)%1 == 0:
        nameall += 'run'+str(i+1)+'.png '

os.system("convert -delay 5 -loop 5 " + nameall +" " + oall)

os.system('rm -rf *.png')
os.system('rm -rf *.h5')
os.system('rm -rf *.yml')
