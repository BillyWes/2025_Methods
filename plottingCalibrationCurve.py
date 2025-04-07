import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

csfont = {'fontname':'Times New Roman',
          'size'    :14}
titlefont = {'fontname':'Times New Roman',
          'size'    :16}

concentrations = [9.5*(10**(-3)), 7.6*(10**(-3)), 5.7*(10**(-3)), 3.8*(10**(-3)), 1.9*(10**(-3)), 9.5*(10**(-4)), 9.5*(10**(-5))]
concentrations = np.array(concentrations)*10**6/253
absorbance448 = [1.814843655, 1.48782897, 1.076797485, 0.723665953, 0.384079099, 0.17765595, 0.025996104]


fit = np.polyfit(concentrations, absorbance448, 1, full=True)
coefficients = fit[0]
m,c = coefficients[0],coefficients[1]
residuals = fit[1]

fit_values = np.polyval([m, c], concentrations)

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams["mathtext.default"] = 'regular'
mpl.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['font.family'] = 'times'

plt.figure(figsize=(8, 6))
plt.grid(True)
plt.scatter(concentrations, absorbance448, color='k')
plt.plot(concentrations, fit_values, color='gray', linestyle='--', label=f"Linear Fit : $ε_{{448 nm}}$ = {int(float('%.3g' % (m*10**6)))} $ \pm$ ${int(float('%.2g' % (residuals[0]*10**6)))}$ $M^{{-1}} cm^{{-1}}$")

plt.xlabel("Concentration (µM)", **titlefont)
plt.xlim(0, 40)
plt.xticks([0, 5, 10, 15, 20, 25, 30, 35, 40],**csfont)
plt.ylabel("Absorbance (448 nm)", **titlefont)
plt.ylim(0, 1.9)
plt.yticks([0.00, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75],**csfont)
plt.title("Merocyanine Aldehyde 5: Concentration v. Absorbance", **titlefont)
plt.legend(prop = {'size': 14}, loc = 'lower right')

plt.savefig("Calibration_Curves/250221_MCy.pdf")
print(residuals)
