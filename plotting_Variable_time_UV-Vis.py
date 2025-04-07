import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
from scipy.optimize import curve_fit
import matplotlib as mpl

csfont = {'fontname':'Times New Roman'}

outPath = 'plottedUV-Vis/'

titles = False
times = np.linspace(0, 480, 145)
timesRange = [0, 480]
wavelengthRange = [199,800]
plotWavelengthRange = [350, 650]
IVdatas = []
MCydatas = []
trackWavelengths = [420, 560]
exclude_times = []
normWL = 280
stoppingTime = 250 
format = 'pdf'

wavelengths = list(range(wavelengthRange[0], wavelengthRange[1]))
wavelengths.reverse()
filename = "UV-Vis/250225Kinetic_Runs_Sccc7-IV_SWAP_ASCII_FORINPUT.csv"  # Change to your actual filename
df = pd.read_csv(filename, skiprows=1)

wave_Abs_alternator = 'wave'
IV_MCy_alternator = 'IV'
excludeCounter = 0

for col in df:
    if excludeCounter > 3:
        if wave_Abs_alternator == 'Abs':
            if IV_MCy_alternator == 'IV':
                IVdatas.append(df[col])
                IV_MCy_alternator = 'MCy'
            else:
                MCydatas.append(df[col])
                IV_MCy_alternator = 'IV'
            wave_Abs_alternator = 'wave'
        else:
            wave_Abs_alternator = 'Abs'
    excludeCounter = excludeCounter+1

def get2Deriv(data):
    dx = 1
    oneD = np.gradient(data, dx)
    return(np.gradient(oneD, dx))

def plot2Deriv():
    tZero = get2Deriv(MCydatas[0])
    tEnd = get2Deriv(MCydatas[len(MCydatas)-2])

    print([wavelengths, tZero, tEnd])

    plt.figure(figsize=(5, 3))

    mpl.rcParams['mathtext.fontset'] = 'stix'
    mpl.rcParams["mathtext.default"] = 'regular'
    mpl.rcParams['font.family'] = 'STIXGeneral'
    plt.rcParams['font.family'] = 'times'

    plt.plot(wavelengths, tZero)
    plt.plot(wavelengths, tEnd)

    plt.show()

def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0))) + b
    return (y)

def plotCurves(dataPack, saveName=''):

    colors = [mcolors.to_hex(c) for c in plt.cm.plasma(np.linspace(0, 1, 1200))]

    plt.figure(figsize=(5, 3))

    downToOneNormFactor = 0

    for i, graphSet in enumerate(dataPack):
        if times[i] not in exclude_times:
            normFactor = (graphSet[(max(wavelengthRange)-min(wavelengthRange)-(normWL-min(wavelengthRange)))])
            normData = (graphSet/normFactor)
            if max(normData[(min(plotWavelengthRange)-min(wavelengthRange)):(max(wavelengthRange)-min(plotWavelengthRange))]) > downToOneNormFactor:
                downToOneNormFactor = max(normData[(min(plotWavelengthRange)-min(wavelengthRange)):(max(wavelengthRange)-min(plotWavelengthRange))])
            
            if times[i] == max(times):
                break

    for i, graphSet in enumerate(dataPack):
        if times[i] not in exclude_times:
            normFactor = (graphSet[(max(wavelengthRange)-min(wavelengthRange)-(normWL-min(wavelengthRange)))])
            if times[i]<=stoppingTime:
                plt.plot(wavelengths, ((graphSet/normFactor)/downToOneNormFactor), color = colors[int(round(1000*(times[i]-min(timesRange))/(stoppingTime-min(timesRange))))], label=times[i], linewidth=0.75)
                print(str(times[i])+' plotted with colour: ' + colors[int(round(1000*(times[i]-min(timesRange))/(stoppingTime-min(timesRange))))] + ' and normFactor: ' + str(normFactor))

            if times[i] == max(times):
                break

    plt.xlabel("Wavelength (nm)", **csfont)
    plt.xlim(350, 650)
    plt.xticks([350,400,450,500,550,600,650], **csfont)
    plt.ylabel("Normalised Absorbance", **csfont)
    plt.ylim(-0.05, 1.15)
    plt.yticks([0,0.2,0.4,0.6,0.8,1], **csfont)
    plt.subplots_adjust(bottom=0.15)
    if titles:
        plt.title("Absorbance Spectra with Varying Time", **csfont)
    plt.grid(False)
    
    if saveName == '':
        plt.show
    else:
        plt.savefig(saveName+'.'+format)
    
def plotSummaryCurves(dataPack, saveName='', quadFit=False, sigFit = False):
    plt.figure(figsize=(5, 3))
    for trackWavelength in trackWavelengths:
        ydata = []
        for i, graphset in enumerate(dataPack):
            normFactor = (graphset[(max(wavelengthRange)-min(wavelengthRange)-(normWL-min(wavelengthRange)))])
            ydata.append(graphset[(max(wavelengthRange)-min(wavelengthRange)-(trackWavelength-min(wavelengthRange)))]/normFactor)
            if times[i] == max(times):
                break
        ydata = ydata/max(ydata)

        #plt.plot(pHs, ydata, 'o', label=str(trackWavelength)+' nm')

        if trackWavelength == 560:
            if quadFit:
                p = np.polyfit(times, ydata, deg=3)
                dp2_d2pH = [p[0]*3*2,p[1]*2]
                pKa = (0-dp2_d2pH[1])/dp2_d2pH[0]

                fitX = np.linspace(7.2,10.8,50)
                fitY = p[0]*fitX**3+p[1]*fitX**2+p[2]*fitX+p[3]

                p9 = np.polyfit(times, ydata, deg=9)

                fit9X = np.linspace(3.5,11,50)
                fit9Y = p9[0]*fitX**9 + p9[1]*fitX**8 + p9[2]*fitX**7 + p9[3]*fitX**6 + p9[4]*fitX**5 + p9[5]*fitX**4 + p9[6]*fitX**3 + p9[7]*fitX**2 + p9[8]*fitX + p9[9]
                plt.plot(fitX, fitY, '-', color = 'gray', label = 'Fit ($pK_a$ = '+str(pKa)+')')
            if sigFit:
                p0 = [max(ydata)-min(ydata), np.median(times),9,min(ydata)]

                fitX = np.linspace(3.5,11,100)

                popt, pcov = curve_fit(sigmoid, times, ydata,p0, method='dogbox')
                fitY = sigmoid(fitX, popt[0], popt[1], popt[2], popt[3])
                
                plt.plot(fitX, fitY, '-', color = 'gray', label = 'Fit ($pK_a$ = '+str(popt[1])+')')
                
            plt.plot(times, ydata, ',', color = 'k')

    plt.xlabel("Time / mins", **csfont)
    plt.xlim(0, 480)
    plt.xticks([0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500],**csfont)
    plt.ylabel("Normalised Absorbance (560 nm)", **csfont)
    plt.ylim(-0.05, 1.15)
    plt.yticks([0,0.2,0.4,0.6,0.8,1], **csfont)
    plt.subplots_adjust(bottom=0.15)
    if titles:
        plt.title("Sccc7-IV + 1 : Absorbance at 560 nm v. Time", **csfont)
    plt.grid(False)

    if saveName == '':
        plt.show
    else:
        plt.savefig(saveName+'.'+format)

plotSummaryCurves(IVdatas, saveName='plottedUV-Vis/NoTitles/KineticsCurve_Sccc7-IV')
plotSummaryCurves(MCydatas, saveName='plottedUV-Vis/NoTitles/KineticsCurve_Sccc7-MCy')
plotCurves(IVdatas, saveName='plottedUV-Vis/NoTitles/Kinetics_Sccc7-IV_NORM')
plotCurves(MCydatas, saveName='plottedUV-Vis/NoTitles/Kinetics_Sccc7-MCy_NORM')

# plot2Deriv()