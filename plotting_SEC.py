import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

# Load the data
file_path = "Analytical HPLC/20241017_c18A CC-Tet-star-+-w-oW-PurifiedB (01).csv" 
blank_path = "Analytical HPLC/20241017_c18A BLANK (01).csv" 
output_path = 'purityPlots/'
name = 'CC-Tet*'
Norm = True
fontsizer = 18
timeLimits = [0, 20]

NormExt = ''
data = pd.read_csv(file_path)
blankData = pd.read_csv(blank_path)

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams["mathtext.default"] = 'regular'
mpl.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['font.family'] = 'times'

if Norm:
    data["Intensity"] = (data["Intensity"])/(max(data["Intensity"]))
else:
    data["Intensity"] = data["Intensity"]
plt.figure(figsize=(5, 5))
plt.plot(data["min"], data["Intensity"], 'k-')

if Norm:
    plt.ylabel("Abs$_{{280}}$", fontsize=fontsizer)
    #plt.ylim(min(data["Intensity"]), max(data["Intensity"])+10/100*max(data["Intensity"]))
    #plt.xlim(min(timeLimits), max(timeLimits))
else:
    plt.ylabel("Abs$_{{280}}$", fontsize=fontsizer)
    #plt.ylim(min(data["Intensity"]), max(data["Intensity"])+10/100*max(data["Intensity"]))
    #plt.xlim(min(timeLimits), max(timeLimits))

plt.xlabel("Time / mins", fontsize=fontsizer)
plt.yticks([])
plt.xlim(min(timeLimits),max(timeLimits))
plt.subplots_adjust(bottom=0.15)
plt.tick_params(axis='both', which='major', labelsize=fontsizer)

if Norm:
    NormExt = 'Norm'
plt.savefig(output_path + name + '_' + NormExt + '.pdf', transparent=True)
