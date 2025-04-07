import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import matplotlib as mpl
from scipy.optimize import curve_fit

#READ: comments to call plotting (either plot all curves or plot base titration curve) at bottom of file

csfont = {'fontname':'Times New Roman'}

outPath = 'plottedUV-Vis/'

titles = False
pHs = [4, 4.5, 5, 5.5, 6, 7.4, 8, 9, 10.145794559, 10.278507906, 10.380016745, 10.531341634, 10.770772905, 11.03738368, 11.41328174, 11.88023563]
pHsRange = [3.5,12.5]
wavelengthRange = [200,701]
plotWavelengthRange = [350, 650]
IVdatas = []
MCydatas = []
trackWavelengths = [420, 560]
exclude_pHs = [4, 4.5, 9.544, 9.176]
normWL = 280
format = 'pdf'

fitExcludeFront = 5
fitExcludeBack = 0

wavelengths = list(range(wavelengthRange[0], wavelengthRange[1]))
wavelengths.reverse()
filename = "UV-Vis/CollectedMCyRunsASCII.csv"  # Change to your actual filename
df = pd.read_csv(filename, skiprows=1)

wave_Abs_alternator = 'wave'
IV_MCy_alternator = 'IV'

for col in df:
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

def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0))) + b
    return (y)

def plotCurves(dataPack, saveName='', zeroSet = False, divisor = 1):

    colors = [mcolors.to_hex(c) for c in plt.cm.coolwarm(np.linspace(0, 1, 1000))]

    plt.figure(figsize=(5, 3))

    downToOneNormFactor = 0
    downToOneTakeOff = 0

    for i, graphSet in enumerate(dataPack):
        if pHs[i] not in exclude_pHs:
            if zeroSet:
                downToOneTakeOff = graphSet[max(wavelengthRange)-max(plotWavelengthRange)]
            normFactor = ((graphSet-downToOneTakeOff)[(max(wavelengthRange)-min(wavelengthRange)-(normWL-min(wavelengthRange)))])
            normData = ((graphSet-downToOneTakeOff)/normFactor)
            if max(normData[(min(plotWavelengthRange)-min(wavelengthRange)):(max(wavelengthRange)-min(plotWavelengthRange))]) > downToOneNormFactor:
                downToOneNormFactor = max(normData[(min(plotWavelengthRange)-min(wavelengthRange)):(max(wavelengthRange)-min(plotWavelengthRange))])
            if pHs[i] == max(pHs):
                break

    for i, graphSet in enumerate(dataPack):
        if pHs[i] not in exclude_pHs:
            if zeroSet:
                downToOneTakeOff = graphSet[max(wavelengthRange)-max(plotWavelengthRange)]
            normFactor = (graphSet[(max(wavelengthRange)-min(wavelengthRange)-(normWL-min(wavelengthRange)))])
            plt.plot(wavelengths, (((graphSet-downToOneTakeOff)/normFactor)/downToOneNormFactor)/divisor, color = colors[int(round(1000*(pHs[i]-min(pHsRange))/(max(pHsRange)-min(pHsRange))))], label=pHs[i])
            print(str(pHs[i])+' plotted with colour: ' + colors[int(round(1000*pHs[i]/14))] + ', normFactor: ' + str(normFactor) + ' and downToOneNormFactor: ' + str(downToOneNormFactor))

            if pHs[i] == max(pHs):
                break

    plt.xlabel("Wavelength (nm)", **csfont)
    plt.xlim(min(plotWavelengthRange), max(plotWavelengthRange))
    plt.xticks([350,400,450,500,550,600,650], **csfont)
    plt.ylabel("Normalised Absorbance", **csfont)
    plt.ylim(-0.05, 1.15)
    plt.yticks([0,0.2,0.4,0.6,0.8,1], **csfont)
    plt.subplots_adjust(bottom=0.15)
    if titles:
        plt.title("Absorbance Spectra with Varying pH", **csfont)
    plt.grid(False)
    
    if saveName == '':
        plt.show
    else:
        plt.savefig(saveName+'.'+format)

def plotSummaryCurves(dataPack, saveName='', quadFit=False, sigFit = False):
    plt.figure(figsize=(5, 3))
    for trackWavelength in trackWavelengths:
        ydata = []
        for graphset in dataPack:
            normFactor = (graphset[(max(wavelengthRange)-min(wavelengthRange)-(normWL-min(wavelengthRange)))])
            ydata.append(graphset[(max(wavelengthRange)-min(wavelengthRange)-(trackWavelength-min(wavelengthRange)))]/normFactor)
            if pHs == max(pHs):
                break
        ydata = ydata/max(ydata)

        #plt.plot(pHs, ydata, 'o', label=str(trackWavelength)+' nm')

        if trackWavelength == 560:
            if quadFit:
                fitting_pHs = (pHs[fitExcludeFront:])
                fitting_ydata = (ydata[fitExcludeFront:])

                p = np.polyfit(fitting_pHs, fitting_ydata, deg=3)
                dp2_d2pH = [p[0]*3*2,p[1]*2]
                pKa = (0-dp2_d2pH[1])/dp2_d2pH[0]

                fitX = np.linspace(7.2,10.8,50)
                fitY = p[0]*fitX**3+p[1]*fitX**2+p[2]*fitX+p[3]

                p9 = np.polyfit(fitting_pHs, fitting_ydata, deg=9)

                fit9X = np.linspace(3.5,11,50)
                fit9Y = p9[0]*fitX**9 + p9[1]*fitX**8 + p9[2]*fitX**7 + p9[3]*fitX**6 + p9[4]*fitX**5 + p9[5]*fitX**4 + p9[6]*fitX**3 + p9[7]*fitX**2 + p9[8]*fitX + p9[9]
                plt.plot(fitX, fitY, '-', color = 'gray', label = 'Fit ($pK_a$ = '+str(pKa)+')')
            if sigFit:
                p0 = [max(ydata)-min(ydata), np.median(pHs),9,min(ydata)]

                fitX = np.linspace(min(pHsRange), max(pHsRange),100)

                popt, pcov = curve_fit(sigmoid, pHs, ydata,p0, method='dogbox')
                fitY = sigmoid(fitX, popt[0], popt[1], popt[2], popt[3])
                
                plt.plot(fitX, fitY, '-', color = 'gray', label = 'Fit ($pK_a$ = '+str(round(popt[1], 2))+')')
                
            plt.plot(pHs, ydata, 'o', color = 'k')

    mpl.rcParams['mathtext.fontset'] = 'stix'
    mpl.rcParams["mathtext.default"] = 'regular'
    mpl.rcParams['font.family'] = 'STIXGeneral'
    plt.rcParams['font.family'] = 'times'

    plt.xlabel("pH", **csfont)
    plt.xlim(min(pHsRange), max(pHsRange))
    plt.xticks([4,5,6,7,8,9,10,11,12],**csfont)
    plt.ylabel("Normalised Absorbance (560 nm)", **csfont)
    plt.ylim(-0.05, 1.15)
    plt.yticks([0,0.2,0.4,0.6,0.8,1], **csfont)
    plt.subplots_adjust(bottom=0.15)
    if titles:
        plt.title("Base Titration of Sccc7-MCy(EKW) + 1 : Absorbance at 560 nm", **csfont)
    plt.grid(False)
    plt.legend(loc=3)

    if saveName == '':
        plt.show
    else:
        plt.savefig(saveName+'.'+format)

# plotCurves(IVdatas, saveName=outPath+'NoTitles/BaseTitration_Sccc7-IV_NORM', zeroSet=True, divisor = 0.903814)
# plotCurves(MCydatas, saveName=outPath+'NoTitles/BaseTitration_Sccc7-MCy_NORM', divisor = 1.05671)

# plotSummaryCurves(IVdatas, saveName=outPath+'NoTitles/BaseTitrationCurve_Sccc7-MCy')
plotSummaryCurves(MCydatas, sigFit = True, saveName=outPath+'Test')


