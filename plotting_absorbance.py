import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random

paths = ['Analytical HPLC/20241017_c18A CC-Tet-star-+-w-oW-PurifiedA (01).csv',
         'Analytical HPLC/20241017_c18A CC-Tet-star-+-w-oW-PurifiedA (02).csv',
         'Analytical HPLC/20241017_c18A CC-Tet-star-+-w-oW-PurifiedB (01).csv',
         'Analytical HPLC/20241017_c18A CC-Tet-star-+-w-oW-PurifiedB (02).csv']         #declare paths
plotTitles = ['Fraction A HPLC Chromatogram of Absorbance', 'Fraction B HPLC Chromatogram of Absorbance']    #declare plot titles: should equal number of paths

normalised = True                   #choose whether graph will be normalised based on highest value
combine = 1                         #choose whether to combine the first n number of graphs onto the same axis (if not required, set to 1)
legend = ['220 nm', '280 nm']       #legend titles for individual lines (when combining graphs onto the same axis)
choice = 'HPLC2'                    #select from preconfigured options - if used, no need to declare the following

xlabel = 'Time / mins'              #declare x-axis title
ylabel = 'Intensity'                #declare y-axis title
expColTitles = ['min', 'Intensity'] #set up expected column titles in the csv

class graphObj:
    def __init__(self, paths, plotTitles, xlabel, ylabel, expColTitles, dimension, normalised, legend):
        self.plotTitles = plotTitles
        self.pathStrings = paths
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.expColTitles = expColTitles
        self.dimension = dimension
        self.normalised = normalised
        self.legend = legend

    def read(self):
        self.timeMins_arr = []
        self.absorbance_arr = []
        normalFactor = 1

        for pathString in self.pathStrings:
            try:
                df = pd.read_csv(pathString, usecols = self.expColTitles)
                if self.normalised:
                    normalFactor = (1/max(df[self.expColTitles[1]].values))
                self.timeMins_arr.append(df[self.expColTitles[0]].values)
                self.absorbance_arr.append(np.multiply((df[self.expColTitles[1]].values), normalFactor))
            except Exception as exc:
                print(f"Missing file path: {exc}")
                return None
    
    def plot(self):
        subCount = 1
        if combine == 1:
            for timeMins in self.timeMins_arr:
                plt.subplot(int(self.dimension[0]),self.dimension[1],subCount)
                plt.plot(timeMins, self.absorbance_arr[int(subCount-1)])
                plt.title(self.plotTitles[int(subCount-1)])
                plt.xlabel(self.xlabel)
                plt.ylabel(self.ylabel)
                subCount = subCount+1
            plt.tight_layout()
            plt.show()
            plt.legend()
        else:
            lineIntMark = 0
            titleIntMark = 0
            for i in range(self.dimension[0]):
                plt.subplot(self.dimension[1], int(self.dimension[0]/combine), i+1)
                lineColour = 'r'
                for i in range(combine):
                    plt.plot(self.timeMins_arr[(subCount-1)], self.absorbance_arr[(subCount-1)], color=lineColour, label=self.legend[lineIntMark])
                    if lineColour == 'r':
                        lineColour='g'
                    elif lineColour == 'g':
                        lineColour='b'
                    elif lineColour == 'b':
                        lineColour = 'c'
                    elif lineColour == 'c':
                        lineColour = 'm'
                    elif lineColour == 'm':
                        lineColour = 'y'
                    elif lineColour == 'y':
                        colourString = '#'
                        for i in 3:
                            colourString = colourString + (hex(random.randint(0,200)))
                        lineColour = colourString
                    subCount = subCount+1
                    if lineIntMark == 0:
                        lineIntMark = 1
                    else:
                        lineIntMark = 0
                plt.title(self.plotTitles[titleIntMark])
                titleIntMark = titleIntMark+1
                plt.xlabel(self.xlabel)
                plt.ylabel(self.ylabel)
            plt.tight_layout()
            plt.legend()
            plt.show()  

def HPLC(paths, plotTitles, dimension, normalised, legend):
    HPLC_plotObj = graphObj(paths, plotTitles, 'Time / mins', 'Intensity', ['min', 'Intensity'], dimension, normalised, legend)
    HPLC_plotObj.read()
    HPLC_plotObj.plot()
    return 'True'

if choice == 'HPLC2':
    print(HPLC(paths, plotTitles, [2,1], normalised, legend))
elif choice == 'HPLC4':
    print(HPLC(paths, plotTitles, [2,2], normalised, legend))
elif choice == '':
    plotObj = graphObj(paths, plotTitles, xlabel, ylabel, expColTitles, [1,(len(paths))], normalised, legend)
    plotObj.read()
    plotObj.plot()