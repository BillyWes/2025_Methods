import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib as mpl
#----------------------------------------------------------------------------------------------------------------
# CC-Tet*_Cy5
startingPeptideConc = 450    #/uM
dilutionFactor = 2
wellNo = 20                 #number of unique wells
fluoWL = '654'              #fluorescence wavelength of interest/nm (654 nm for Cy5, 591 nm for MCy5)
fluoCol = 'L'
oligoState = 4              #peptide oligomeric state
fluoPath = 'Plate_Reader/WW-Cy5_2uM-40-CC-Tet-star-Cy5-Fluo.csv'
FPPath = 'Plate_Reader/WW-Cy5_2uM-40-CC-Tet-star-Cy5-FP.csv'
savePath = 'plotted-plates/Figs/CC-Tet*/CC-Tet*_Cy5_'
excludeFluoIndeces = [4]    #exclude indeces of final arrays by index - DESCENDING order!
excludeFPIndeces = [10, 3]
fluorFitted = False
FPFitted = False
fontsizer = 16

fluoYRange = []
fluoNormYRange = []
FPYRange = []
FPNormYRange = []
#----------------------------------------------------------------------------------------------------------------
# sc-CC-7-IV_Cy5
# startingPeptideConc = 20    #/uM
# dilutionFactor = 2
# wellNo = 20                 #number of unique wells
# fluoWL = '654'              #fluorescence wavelength of interest/nm (654 nm for Cy5, 591 nm for MCy5)
# fluoCol = 'L'
# oligoState = 1              #peptide oligomeric state
# fluoPath = 'Plate_Reader/Billy-250214_Cy5CH3_2uM_Sccc7-IV_Gain-2000.CSV'
# FPPath = 'Plate_Reader/Billy-250214_Cy5CH3_2uM_Sccc7-IV_FP.CSV'
# savePath = 'plotted-plates/Figs/sc-CC-7-IV/Sccc7-IV_Cy5_'
# excludeFluoIndeces = []    #exclude indeces of final arrays by index - DESCENDING order!
# excludeFPIndeces = []
# fluorFitted = False
# FPFitted = False
# fontsizer = 16
#----------------------------------------------------------------------------------------------------------------
# # sc-CC-7-IV_MCy5
# startingPeptideConc = 20    #/uM
# dilutionFactor = 2
# wellNo = 20                 #number of unique wells
# fluoWL = '591'              #fluorescence wavelength of interest/nm (654 nm for Cy5, 591 nm for MCy5)
# fluoCol = 'L'
# oligoState = 1              #peptide oligomeric state
# fluoPath = 'Plate_Reader/250226_MeroCy5_Sccc7-IV_binding.CSV'
# FPPath = ''
# savePath = 'plotted-plates/Figs/sc-CC-7-IV/Sccc7-IV_MCy5_'
# excludeFluoIndeces = [9, 5]    #exclude indeces of final arrays by index - DESCENDING order!
# excludeFPIndeces = []
# fluorFitted = False
# FPFitted = False
# fontsizer = 16
#----------------------------------------------------------------------------------------------------------------
# # sc-CC-7-MCy(EKW)_MCy5
# startingPeptideConc = 20    #/uM
# dilutionFactor = 2
# wellNo = 20                 #number of unique wells
# fluoWL = '591'              #fluorescence wavelength of interest/nm (654 nm for Cy5, 591 nm for MCy5)
# fluoCol = 'L'
# oligoState = 1              #peptide oligomeric state
# fluoPath = 'Plate_Reader/250226_MeroCy5_Sccc7-MCy_binding.CSV'
# FPPath = ''
# savePath = 'plotted-plates/Figs/sc-CC-7-MCy(EKW)/Sccc7-MCy(EKW)_MCy5_'
# excludeFluoIndeces = [16]    #exclude indeces of final arrays by index - DESCENDING order!
# excludeFPIndeces = []
# fluorFitted = False
# FPFitted = False
# fontsizer = 16
#----------------------------------------------------------------------------------------------------------------
# # sc-CC-7-MCy(EKW)_Cy5
# startingPeptideConc = 20    #/uM
# dilutionFactor = 2
# wellNo = 20                 #number of unique wells
# fluoWL = '654'              #fluorescence wavelength of interest/nm (654 nm for Cy5, 591 nm for MCy5)
# fluoCol = 'L'
# oligoState = 1              #peptide oligomeric state
# fluoPath = 'Plate_Reader/Billy-250212_Cy5CH3_2uM_Sccc7-MCy(EKW)_Gain-2000.CSV'
# FPPath = 'Plate_Reader/Billy-250212_Cy5CH3_2uM_Sccc7-MCy(EKW)_FP.CSV'
# savePath = 'plotted-plates/Figs/sc-CC-7-MCy(EKW)/Sccc7-MCy(EKW)_Cy5_'
# excludeFluoIndeces = [16, 9]    #exclude indeces of final arrays by index - DESCENDING order!
# excludeFPIndeces = []
# fluorFitted = False
# FPFitted = False
# fontsizer = 16
#----------------------------------------------------------------------------------------------------------------

def quadratic_binding(x, Bmax, Kd):
    return (Bmax * x) / (Kd + x)

def fit_binding_curve(x_data, y_data, print_ = False):
    popt, _ = curve_fit(quadratic_binding, x_data, y_data, p0=[max(y_data), np.median(x_data)])
    if print_:
        return _
    else:
        return popt

def getPeptideConc(choice):
    wellPeptideConcs = []
    if choice == 'barrel':
        newConc = (startingPeptideConc)/oligoState
    else:
        newConc = (startingPeptideConc)
    for i in range(wellNo-1):
     wellPeptideConcs.append(newConc)
     newConc = newConc/dilutionFactor
    wellPeptideConcs.append(0)
    return(wellPeptideConcs)

def getFluoAbsor(choice):
    f = pd.read_csv(fluoPath, skiprows = 6)
    columnInterest = f[fluoWL]
    means = []
    errors = []
    for i in range(wellNo):
       o=20-(i+1)
       means.append((columnInterest[o]+columnInterest[o+20]+columnInterest[o+40]+columnInterest[o+60])/4)
       errors.append(np.nanstd([columnInterest[o], columnInterest[o+20], columnInterest[o+40], columnInterest[o+60]], axis=0))
    if choice == 'means':
        return(means)
    else:
        return(errors)   
    
def getFPs(choice):
    f = pd.read_csv(FPPath, skiprows = 5, header = 0)
    columnInterest = f[' Polarization based on Raw Data (F: 590-50 / F: 675-50)']
    means = []
    errors = []
    for i in range(wellNo):
        o=20-(i+1)
        means.append((columnInterest[o]+columnInterest[o+20]+columnInterest[o+40]+columnInterest[o+60])/4)
        errors.append(np.nanstd([columnInterest[o], columnInterest[o+20], columnInterest[o+40], columnInterest[o+60]], axis=0))
    if choice == 'means':
        return(means)
    else:
        return(errors)


def plotFluoInt(deRepeater=True):

    stdFigAx = ''

    plt.figure(figsize=(5, 5))

    mpl.rcParams['mathtext.fontset'] = 'stix'
    mpl.rcParams["mathtext.default"] = 'regular'
    mpl.rcParams['font.family'] = 'STIXGeneral'
    plt.rcParams['font.family'] = 'times'

    barrelConc = np.asarray(getPeptideConc('barrel'))
    fluoAbsor = np.asarray(getFluoAbsor('means'))
    errors = np.asarray(getFluoAbsor(''))
    if deRepeater:
        normStat = '_norm'
    else:
        normStat = '_raw'

    for excludeIndex in excludeFluoIndeces:
            barrelConc = np.delete(barrelConc, excludeIndex)
            fluoAbsor = np.delete(fluoAbsor, excludeIndex)
            errors = np.delete(errors, excludeIndex)

    if deRepeater:
        errors = errors / max(fluoAbsor)
        fluoAbsor = fluoAbsor / max(fluoAbsor)

        if fluoNormYRange != []:
            plt.ylim(min(fluoNormYRange), max(fluoNormYRange))
            stdFigAx = '_STD'

    else:
        fluoAbsor = fluoAbsor/1000
        errors = errors/1000

        if fluoYRange != []:
            plt.ylim(min(fluoYRange), max(fluoYRange))
            stdFigAx = '_STD'

    if fluorFitted:
        Bmax, Kd = fit_binding_curve(barrelConc, fluoAbsor)
        errorsOf=fit_binding_curve(barrelConc, fluoAbsor, print_=True)
        errorOf=np.sqrt(np.diag(errorsOf))
        fit_x = np.linspace(min(barrelConc), max(barrelConc), 100)
        fit_y = quadratic_binding(fit_x, Bmax, Kd)
        plt.plot(fit_x, fit_y, label=f'Fit (K$_d$ = {Kd:.2f} $\\pm$ {errorOf[0]:.2f} µM)', color = 'gray')
        plt.legend(fontsize=fontsizer)

    plt.errorbar(barrelConc, fluoAbsor, errors, fmt='.', elinewidth=1, ecolor='k', color='k')
    if oligoState == 1:
        plt.xlabel('Protein Concentration / µM', fontsize=fontsizer)
        plt.xticks([0,5,10,15,20])
    else:
        plt.xlabel('Barrel Concentration / µM', fontsize=fontsizer)
    plt.xticks(fontsize=fontsizer)
    plt.yticks(fontsize=fontsizer)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    if fluoAbsor.max() > 1.2:
        plt.ylabel('FI$_{}$'.format(fluoWL[0])+'$_{}$'.format(fluoWL[1])+'$_{}$'.format(fluoWL[2]) + ' / 10$^{3}$', fontsize=fontsizer)
    else:
        plt.ylabel('Normalised FI$_{}$'.format(fluoWL[0])+'$_{}$'.format(fluoWL[1])+'$_{}$'.format(fluoWL[2]), fontsize=fontsizer)
    if fluorFitted:
        plt.legend(loc=4,fontsize=fontsizer)

    #plt.show()
    plt.savefig(savePath+'Fluorescence' + normStat + stdFigAx + '.pdf', transparent=True)
    plt.clf()

    if deRepeater:
        plotFluoInt(False)

def plotFP(deRepeater=True):

    plt.figure(figsize=(5, 5))

    mpl.rcParams['mathtext.fontset'] = 'stix'
    mpl.rcParams["mathtext.default"] = 'regular'
    mpl.rcParams['font.family'] = 'STIXGeneral'
    plt.rcParams['font.family'] = 'times'

    barrelConc = np.asarray(getPeptideConc('barrel'))
    FPs = np.asarray(getFPs('means'))
    errors = np.asarray(getFPs(''))

    if deRepeater:
        normStat = '_norm'
    else:
        normStat = '_raw'

    for excludeIndex in excludeFPIndeces:
        barrelConc = np.delete(barrelConc, excludeIndex)
        FPs = np.delete(FPs, excludeIndex)
        errors = np.delete(errors, excludeIndex)

    if FPYRange != []:
        plt.ylim(min(FPYRange), max(FPYRange))
        stdFigAx = '_STD'

    if deRepeater:
        errors = errors / max(FPs)
        FPs = FPs / max(FPs)

        if FPNormYRange != []:
            plt.ylim(min(FPNormYRange), max(FPNormYRange))
            stdFigAx = '_STD'

    plt.errorbar(barrelConc, FPs, errors, fmt='.', elinewidth=1, ecolor='k', color='k')

    if FPFitted:
        Bmax, Kd = fit_binding_curve(barrelConc, FPs)
        errorsOf=fit_binding_curve(barrelConc, FPs, print_=True)
        errorOf=np.sqrt(np.diag(errorsOf))
        fit_x = np.linspace(min(barrelConc), max(barrelConc), 100)
        fit_y = quadratic_binding(fit_x, Bmax, Kd)
        plt.plot(fit_x, fit_y, label=f'Fit (K$_d$ = {Kd:.2f} $ \\pm$ {errorOf[0]:.2f} µM)', color = 'gray')
    
    plt.xticks(fontsize=fontsizer)
    plt.yticks(fontsize=fontsizer)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    if oligoState == 1:
        plt.xlabel('Protein Concentration / µM', fontsize=fontsizer)
        plt.xticks([0,5,10,15,20])
    else:
        plt.xlabel('Barrel Concentration / µM', fontsize=fontsizer)
    if FPs.max() > 1.2:
        plt.ylabel('FP', fontsize=fontsizer)
    else:
        plt.ylabel('Normalised FP', fontsize=fontsizer)
    if FPFitted:
       plt.legend(loc=4, fontsize=fontsizer)
    
    #plt.show()
    plt.savefig(savePath + 'FP' + normStat + stdFigAx + '.pdf', transparent=True)
    plt.clf()

    if deRepeater:
        plotFP(False)

def replotAllCharts():

    global startingPeptideConc, dilutionFactor, wellNo, fluoWL, fluoCol, oligoState, fluoPath, FPPath, savePath, excludeFluoIndeces, excludeFPIndeces, fluorFitted, FPFitted, fontsizer, fluoYRange, fluoNormYRange, FPYRange, FPNormYRange

    fluoYRange = [-5, 250]
    fluoNormYRange = [-0.05, 1.2]
    FPYRange = [-50, 200]
    FPNormYRange = [-0.3, 1.2]

    # CC-Tet*_Cy5
    startingPeptideConc = 450    #/uM
    dilutionFactor = 2
    wellNo = 20                 #number of unique wells
    fluoWL = '654'              #fluorescence wavelength of interest/nm (654 nm for Cy5, 591 nm for MCy5)
    fluoCol = 'L'
    oligoState = 4              #peptide oligomeric state
    fluoPath = 'Plate_Reader/WW-Cy5_2uM-40-CC-Tet-star-Cy5-Fluo.csv'
    FPPath = 'Plate_Reader/WW-Cy5_2uM-40-CC-Tet-star-Cy5-FP.csv'
    savePath = 'plotted-plates/Figs_Replots/CC-Tet*/CC-Tet*_Cy5_'
    excludeFluoIndeces = [4]    #exclude indeces of final arrays by index - DESCENDING order!
    excludeFPIndeces = [10, 3]
    fluorFitted = False
    FPFitted = False
    fontsizer = 16

    plotFluoInt()
    plotFP()

    # sc-CC-7-IV_Cy5
    startingPeptideConc = 20    #/uM
    dilutionFactor = 2
    wellNo = 20                 #number of unique wells
    fluoWL = '654'              #fluorescence wavelength of interest/nm (654 nm for Cy5, 591 nm for MCy5)
    fluoCol = 'L'
    oligoState = 1              #peptide oligomeric state
    fluoPath = 'Plate_Reader/Billy-250214_Cy5CH3_2uM_Sccc7-IV_Gain-2000.CSV'
    FPPath = 'Plate_Reader/Billy-250214_Cy5CH3_2uM_Sccc7-IV_FP.CSV'
    savePath = 'plotted-plates/Figs_Replots/sc-CC-7-IV/Sccc7-IV_Cy5_'
    excludeFluoIndeces = []    #exclude indeces of final arrays by index - DESCENDING order!
    excludeFPIndeces = []
    fluorFitted = False
    FPFitted = True
    fontsizer = 16

    plotFluoInt()
    plotFP()

    # sc-CC-7-MCy(EKW)_Cy5
    startingPeptideConc = 20    #/uM
    dilutionFactor = 2
    wellNo = 20                 #number of unique wells
    fluoWL = '654'              #fluorescence wavelength of interest/nm (654 nm for Cy5, 591 nm for MCy5)
    fluoCol = 'L'
    oligoState = 1              #peptide oligomeric state
    fluoPath = 'Plate_Reader/Billy-250212_Cy5CH3_2uM_Sccc7-MCy(EKW)_Gain-2000.CSV'
    FPPath = 'Plate_Reader/Billy-250212_Cy5CH3_2uM_Sccc7-MCy(EKW)_FP.CSV'
    savePath = 'plotted-plates/Figs_Replots/sc-CC-7-MCy(EKW)/Sccc7-MCy(EKW)_Cy5_'
    excludeFluoIndeces = [16, 9]    #exclude indeces of final arrays by index - DESCENDING order!
    excludeFPIndeces = []
    fluorFitted = False
    FPFitted = True
    fontsizer = 16

    plotFluoInt()
    plotFP()

    fluoYRange = [0, 30]
    fluoNormYRange = [0, 1.25]

    # sc-CC-7-IV_MCy5
    startingPeptideConc = 20    #/uM
    dilutionFactor = 2
    wellNo = 20                 #number of unique wells
    fluoWL = '591'              #fluorescence wavelength of interest/nm (654 nm for Cy5, 591 nm for MCy5)
    fluoCol = 'L'
    oligoState = 1              #peptide oligomeric state
    fluoPath = 'Plate_Reader/250226_MeroCy5_Sccc7-IV_binding.CSV'
    FPPath = ''
    savePath = 'plotted-plates/Figs_Replots/sc-CC-7-IV/Sccc7-IV_MCy5_'
    excludeFluoIndeces = [9, 5]    #exclude indeces of final arrays by index - DESCENDING order!
    excludeFPIndeces = []
    fluorFitted = False
    FPFitted = False
    fontsizer = 16

    plotFluoInt()

    # sc-CC-7-MCy(EKW)_MCy5
    startingPeptideConc = 20    #/uM
    dilutionFactor = 2
    wellNo = 20                 #number of unique wells
    fluoWL = '591'              #fluorescence wavelength of interest/nm (654 nm for Cy5, 591 nm for MCy5)
    fluoCol = 'L'
    oligoState = 1              #peptide oligomeric state
    fluoPath = 'Plate_Reader/250226_MeroCy5_Sccc7-MCy_binding.CSV'
    FPPath = ''
    savePath = 'plotted-plates/Figs_Replots/sc-CC-7-MCy(EKW)/Sccc7-MCy(EKW)_MCy5_'
    excludeFluoIndeces = [16]    #exclude indeces of final arrays by index - DESCENDING order!
    excludeFPIndeces = []
    fluorFitted = True
    FPFitted = False
    fontsizer = 16

    plotFluoInt()

replotAllCharts()

