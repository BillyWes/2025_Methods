import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

filename = 'UV-Vis/Solvents/240225_MCy_18-75uM_n-butylamine_MeOH_PBS_7-4_10mM.csv'
df = pd.read_csv(filename, skiprows=1)

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams["mathtext.default"] = 'regular'
mpl.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['font.family'] = 'times'

print(df)

plt.figure(figsize=(5, 3))
maxAbs = max(df['Abs'])

plt.plot(df['Wavelength (nm)'], df['Abs'], 'k-', label = '$\\lambda_{{max}}$ = ' + str(maxAbs) + ' nm')
plt.xlabel("Wavelength / nm")
plt.ylabel("Absorbance")
plt.legend()

plt.show()