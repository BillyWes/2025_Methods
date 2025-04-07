import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.interpolate import make_interp_spline, BSpline

#xydata

wavelengths =   [190, 191, 192.5, 195, 197, 200, 202, 205, 208, 210, 211, 214, 215, 217, 220, 222, 225, 230, 234, 238, 240, 250]

aHelix =        [748, 769, 733, 643, 443, 143, 0, -250, -326, -324, -321, -310, -314, -331, -353, -357, -324, -219, -114, -43, -33, 0]
bStrand =       [224, 253, 300, 319, 300, 243, 193, 57, -47, -108, -121, -164, -179, -184, -157, -138, -114, -64, -36, -14, 7, 0]
randomCoil =    [-322, -347, -375, -410, -419, -364, -256, -145, -34, -14, 0, 35, 41, 46, 44, 39, 27, 8, 0, -1.4, -1.5, 0]

baseline = np.linspace(0,0,22)

def plotSeries(x, y, colorer, namer, linestyler):
    smoothX = np.linspace(min(x), max(x), 100)
    interpSmoother = make_interp_spline(x, y, k=3)
    smoothY = interpSmoother(smoothX)

    plt.plot(smoothX, smoothY/10, color = colorer, linestyle = linestyler, label = namer, linewidth = 1.6)

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams["mathtext.default"] = 'regular'
mpl.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['font.family'] = 'times'

plt.figure(figsize=(5, 6))

plt.plot(wavelengths, baseline, color = 'grey', linestyle = 'dotted', linewidth = 0.7)

plotSeries(wavelengths, aHelix, '#62D5A0', '$\\alpha$-Helix', 'solid')
plotSeries(wavelengths, bStrand, '#CF8AD3', '$\\beta$-Strand', 'dashed')
plotSeries(wavelengths, randomCoil, '#2452C0', 'Random Coil', 'dashdot')

plt.subplots_adjust(bottom=0.15)
plt.legend()
plt.xlabel('Wavelength / nm')
plt.ylabel('$MRE$' + ' ' + '$/$' + ' ' + '$ deg$' + ' ' + '$ cm^2$' + ' ' + '$ dmol^{-1}$' + ' ' + '$ 10^{-3}$')
plt.xlim([190, 250])

plt.savefig('plottedCDs/exemplarCD.png')
plt.savefig('plottedCDs/exemplarCD.pdf')