import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

csfont = {'fontname':'Times New Roman'}

noOfTurns = 10
ResiPerTurn = 3.6
risePerResi = 1.5

helix1y = []
helix1x = []

degreesPerResi = 360/ResiPerTurn

def getXVal(angle):
    localTurnCount = 0
    while angle>360:
        angle = angle-360
        localTurnCount = localTurnCount+1
    if localTurnCount > noOfTurns:
        return(0.1)
    else:
        return(angle)

#generate angles (x values)

i = 0
iHad = 0
while iHad%1 != 0.1:
    iHad=getXVal(i)
    helix1x.append(iHad)
    i=i+degreesPerResi
helix1x = helix1x[:-1]

o = 0
while len(helix1y)<len(helix1x):
    helix1y.append(o)
    o = o+risePerResi

dResiX = helix1x[6::7]
dResiY = helix1y[6::7]

aResiX = helix1x[2::7]
aResiY = helix1y[2::7]

fig = plt.figure()
ax = fig.add_subplot(111)

plt.plot(helix1x, helix1y,'.', color = '#86EEC7')
plt.plot(dResiX, dResiY,'o', color = '#50B78B')
plt.plot(aResiX, aResiY,'o', color = '#50B78B')
plt.xlabel('angle / °', **csfont)
plt.ylabel('helix height / Å', **csfont)

#plt.xticks((0,60,120,180,240,300,360), **csfont)
#plt.yticks((0,10,20,30,40), **csfont)

plt.xticks([])
plt.yticks([])

letterVals = ['c', 'b', 'a', 'g', 'f', 'e', 'd']

for i, angleic in enumerate(helix1x):
    annoteChoice = letterVals[((i)%7)]
    colourChoice = 'gray'
    if annoteChoice == 'a' or 'd':
        colourChoice == 'black'
    #ax.annotate(('  '+annoteChoice), (helix1x[i], helix1y[i]), color = colourChoice, **csfont)

fig.savefig('helicalNet1Unlabelled.png', transparent=True)