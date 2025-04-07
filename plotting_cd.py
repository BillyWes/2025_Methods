import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib as mpl

#----------------------------------------------------------------------------------------------------------------
# parentPath = 'CD/242910_CC-Tet-star/'             #folder containing CD output files
# savePath = 'plottedCDs/242910_CC-Tet-star/'       #REMEMBER to trail with /
# saveName = 'CC-Tet-star_Replot_1'

# blankPath = '242910_blank_0uM_CC-Tet-star-amides_1mm_20C_PBS74.txt'   #name of text file for associated trace
# coolPath =  '242910_cool_100uM_CC-Tet-star-amides_1mm_20C_PBS74.txt'
# meltPath =  '242910_melt_100uM_CC-Tet-star-amides_1mm_20C_PBS74.txt'
# postPath =  '242910_post_100uM_CC-Tet-star-amides_1mm_20C_PBS74.txt'
# prePath =   '242910_pre_100uM_CC-Tet-star-amides_1mm_20C_PBS74.txt'
# fontsizer = 24

# peptideID = 'CC-Tet-star'                  #name of peptide in Peptide_Sheet.csv - if not present, update
# peptideConc = 100                     #concentration of peptide / µM
# pathLength = 1                          #path length / mm
# meltCoolWavelength = 222                #wavelength for which melt/cool moar ellipticity was measured
# wavelengthRange = [200, 250]            #range for displayed wavelengths on the MeltCool chart
# molarElipticityRange = [-45,35]
# molarElipticityRangeMeltCool = [-40, -20]
# percentage = True                       #y-scale percentage helicity (*10^-3)

# lineAt700 = True                        #toggle a line demarkating unreliable CD signals at HT of >700 V
#----------------------------------------------------------------------------------------------------------------
parentPath = 'CD/250206_Sccc7-IV/'             #folder containing CD output files
savePath = 'plottedCDs/250206_Sccc7-IV/'       #REMEMBER to trail with /
saveName = 'Sccc7-IV_Replot'

blankPath = '20250206_blank_0uM_Sccc7-IV_368amides_1mm_20C_10mM_PBS.txt'   #name of text file for associated trace
coolPath =  '20250206_cool_7-2uM_Sccc7-IV_368amides_1mm_20C_10mM_PBS.txt'
meltPath =  '20250206_melt_7-2uM_Sccc7-IV_368amides_1mm_20C_10mM_PBS.txt'
postPath =  '20250206_post_7-2uM_Sccc7-IV_368amides_1mm_20C_10mM_PBS.txt'
prePath =   '20250206_pre_7-2uM_Sccc7-IV_368amides_1mm_20C_10mM_PBS.txt'
fontsizer = 24

peptideID = 'Sccc7-IV'                  #name of peptide in Peptide_Sheet.csv - if not present, update
peptideConc = 4.02                     #concentration of peptide / µM
pathLength = 1                          #path length / mm
meltCoolWavelength = 222                #wavelength for which melt/cool moar ellipticity was measured
wavelengthRange = [200, 250]            #range for displayed wavelengths on the MeltCool chart
molarElipticityRange = [-45,35]
molarElipticityRangeMeltCool = [-40, -20]
percentage = True                       #y-scale percentage helicity (*10^-3)

lineAt700 = True                        #toggle a line demarkating unreliable CD signals at HT of >700 V
#----------------------------------------------------------------------------------------------------------------
# parentPath = 'CD/250213_Sccc7-MCy(EKW)/'             #folder containing CD output files
# savePath = 'plottedCDs/250213_Sccc7-MCy(EKW)/'       #REMEMBER to trail with /
# saveName = 'Sccc7-MCy(EKW)_Replot_1'

# blankPath = '240416_blank_0uM_Sccc7-MCy(EKW)_394amides_1mm_20C_PBS-pH=7-4.txt'   #name of text file for associated trace
# coolPath =  '240416_cool_6-75uM_Sccc7-MCy(EKW)_394amides_1mm_20C_PBS-pH=7-4.txt'
# meltPath =  '240416_melt_6-75uM_Sccc7-MCy(EKW)_394amides_1mm_20C_PBS-pH=7-4.txt'
# postPath =  '240416_post_6-75uM_Sccc7-MCy(EKW)_394amides_1mm_20C_PBS-pH=7-4.txt'
# prePath =   '240416_pre_6-75uM_Sccc7-MCy(EKW)_394amides_1mm_20C_PBS-pH=7-4.txt'
# fontsizer = 24

# peptideID = 'Sccc7-MCy(EKW)'                  #name of peptide in Peptide_Sheet.csv - if not present, update
# peptideConc = 6.75                     #concentration of peptide / µM
# pathLength = 1                          #path length / mm
# meltCoolWavelength = 222                #wavelength for which melt/cool moar ellipticity was measured
# wavelengthRange = [200, 250]            #range for displayed wavelengths on the MeltCool chart
# molarElipticityRange = [-45,35]
# molarElipticityRangeMeltCool = [-40, -20]
# percentage = True                       #y-scale percentage helicity (*10^-3)

# lineAt700 = True                        #toggle a line demarkating unreliable CD signals at HT of >700 V
#----------------------------------------------------------------------------------------------------------------
def readFileBPP(path):                  #takes in a blank, post or pre path and returns a table with labelled headers
    f = pd.read_csv(path, sep="\t", skiprows=19)
    f.columns = ['Wavelength/nm','CD/mdeg','HT/V']
    return(f)

def readFileMeltCool(path):             #takes in a melt or cool path and returns a table with labelled headers
    f = pd.read_csv(path, sep="\t", skiprows=19)
    f.columns = ['Temperature/C','CD/mdeg','HT/V']
    return(f)

def getConstants():                     #used to extract constants from the peptide cheat sheet using the constant 'peptideID'. Returns a table of constants with labelled headers
    f = pd.read_csv('Peptide_Sheet.csv', skiprows = 0)
    f.columns = ['peptideID','molWeight','length/res','TFA/count', 'amBond/count']
    return(f[f['peptideID']==peptideID])

def getMRE(pathChoice, returnChoice, defChoice = True):   #takes in a path choice (e.g. meltPath) and returns an array of values corresponding to the name 'returnChoice'
    blankTable = readFileBPP(parentPath+blankPath)
    Table = readFileMeltCool(parentPath+pathChoice)
    percentageFactor = 1
    if percentage:
        percentageFactor = 1000 
    if returnChoice == 'MCCDs':
        Table = readFileMeltCool(parentPath+pathChoice)
        blankCD = blankTable.iloc[37,1]
        CDs = Table['CD/mdeg']
        MREs = []
        for CD in CDs:
            MREs.append(((CD-blankCD)/(peptideConcM*amBondConst*pathLength)/percentageFactor))
        return(MREs)
    elif returnChoice == 'PPCDs':
        Table = readFileBPP(parentPath+pathChoice)
        CDs = Table['CD/mdeg']
        MREs = []
        indexer = 0
        for CD in CDs:
            MREs.append(((CD-blankTable.iloc[indexer, 1])/(peptideConcM*amBondConst*pathLength))/percentageFactor)
        return(MREs)
    elif returnChoice == 'MCHT/V':
        Table = readFileMeltCool(parentPath+pathChoice)
        return(Table['HT/V'])
    elif returnChoice == 'PPHT/V':
        Table = readFileBPP(parentPath+pathChoice)
        return(Table['HT/V'])
    elif returnChoice == 'Temperature/C':
        Table = readFileMeltCool(parentPath+pathChoice)
        return(Table[returnChoice])
    elif returnChoice == 'Wavelength/nm':
        Table = readFileBPP(parentPath+pathChoice)
        return(Table[returnChoice])

def plotSaveMeltCool():                 #plots and saves a chart for the melting and cooling curves

    figMeltCool = plt.figure(figsize=(8, 8))
    gs = GridSpec(2, 1, height_ratios=[4, 1])
    percentageString = ''
    if percentage:
        percentageString = '$ 10^{-3}$'

    mpl.rcParams['mathtext.fontset'] = 'stix'
    mpl.rcParams["mathtext.default"] = 'regular'
    mpl.rcParams['font.family'] = 'STIXGeneral'
    plt.rcParams['font.family'] = 'times'

    ax1 = plt.subplot(gs[0])

    ax1.plot(getMRE(meltPath, 'Temperature/C'), getMRE(meltPath, 'MCCDs'), 'k-', label = 'Melt')
    ax1.plot(getMRE(coolPath, 'Temperature/C'), getMRE(coolPath, 'MCCDs'), 'k--', label = 'Cool')

    ax1.legend(fontsize=fontsizer, loc=4)
    ax1.set_ylabel('$MRE_{222}$' + ' ' + '$/$' + ' ' + '$ deg$' + ' ' + '$ cm^2$' + ' ' + '$ dmol^{-1}$' + ' ' + '$ res^{-1}$' + ' ' + percentageString, fontsize=fontsizer)
    ax1.set_ybound(molarElipticityRangeMeltCool)
    ax1.set_yticks([-40, -30, -20])
    ax1.tick_params(axis='both', which='major', labelsize=fontsizer)
    #ax1.set_xlabel('$Temperature$' + ' ' + '$ /$' + ' ' + '$ ^oC$', fontsize=fontsizer)

    ax2 = plt.subplot(gs[1])
    ax2.plot(getMRE(meltPath, 'Temperature/C'), getMRE(meltPath, 'MCHT/V'), 'k-', label = 'Melt')
    ax2.plot(getMRE(coolPath, 'Temperature/C'), getMRE(coolPath, 'MCHT/V'), 'k--', label = 'Cool')

    ax2.set_xlabel('$Temperature$' + ' ' + '$ /$' + ' ' + '$ ^oC$', fontsize=fontsizer)
    ax2.set_ylabel('$HT $' + ' ' + '$/ $' + ' ' + '$V$', fontsize=fontsizer)
    ax2.set_ybound(375, 425)
    ax2.set_yticks([375, 425])

    ax2.tick_params(axis='both', which='major', labelsize=fontsizer)

    figMeltCool.tight_layout()

    figMeltCool.savefig(savePath+saveName+'_MeltCool.pdf', transparent=True)

def plotSavePrePost():                  #plots and saves a chart for the pre and post folding curves

    figPrePost = plt.figure(figsize=(8, 8))
    gs = GridSpec(2, 1, height_ratios=[4, 1])
    percentageString = ''
    if percentage:
        percentageString = '$ 10^{-3}$'

    mpl.rcParams['mathtext.fontset'] = 'stix'
    mpl.rcParams["mathtext.default"] = 'regular'
    mpl.rcParams['font.family'] = 'STIXGeneral'
    plt.rcParams['font.family'] = 'times'

    ax1 = plt.subplot(gs[0])

    ax1.plot(getMRE(prePath, 'Wavelength/nm'), getMRE(prePath, 'PPCDs'), 'k-', label = 'Pre')
    ax1.plot(getMRE(postPath, 'Wavelength/nm'), getMRE(postPath, 'PPCDs'), 'k--', label = 'Post')
    ax1.legend(fontsize=fontsizer)
    ax1.set_ylabel('$MRE$' + ' ' + '$/$' + ' ' + '$ deg$' + ' ' + '$ cm^2$' + ' ' + '$ dmol^{-1}$' + ' ' + '$ res^{-1}$' + ' ' + percentageString, fontsize=fontsizer)
    ax1.set_xbound(wavelengthRange)
    ax1.set_ybound(molarElipticityRange)

    ax1.set_yticks([-40, -20, 0, 20])
    ax1.tick_params(axis='both', which='major', labelsize=fontsizer)
    #ax1.set_xlabel('$Wavelength$' + ' ' + '$ /$' + ' ' + '$ nm$', fontsize=fontsizer)

    ax2 = plt.subplot(gs[1])
    if lineAt700:
        ax2.plot(wavelengthRange,[700,700], color='#c9c9c9', linestyle = 'dotted')

    ax2.plot(getMRE(prePath, 'Wavelength/nm'), getMRE(prePath, 'PPHT/V'), 'k-', label = 'Pre')
    ax2.plot(getMRE(postPath, 'Wavelength/nm'), getMRE(postPath, 'PPHT/V'), 'k--', label = 'Post')
    ax2.set_xlabel('$Wavelength$' + ' ' + '$ /$' + ' ' + '$ nm$', fontsize=fontsizer)
    ax2.set_ylabel('$HT $' + ' ' + '$/ $' + ' ' + '$V$', fontsize=fontsizer)
    ax2.set_xbound(wavelengthRange)

    ax2.tick_params(axis='both', which='major', labelsize=fontsizer)

    figPrePost.tight_layout()

    print(getMRE(prePath, 'PPCDs'))
    print('MRE_208 nm : ' + getMRE(prePath, 'PPCDs')[208])
    print('MRE_222 nm : ' + getMRE(prePath, 'PPCDs')[222])

    figPrePost.savefig(savePath+saveName+'_PrePost.pdf', transparent=True)

amBondConst = (getConstants())['amBond/count']
peptideConcM = peptideConc*10**(-6)

plotSaveMeltCool()
plotSavePrePost()