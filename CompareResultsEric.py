import numpy as np
import matplotlib.pyplot as plt

Eric = [[7,8,10,12,15,20,25,30],[7.90E-4,1.17E-4,1.35E-5,4.17E-6,1.55E-6,7.34E-7,5.47E-7,4.85E-7]]

MLE100 = [3.56E-4,2.76E-5,6.26E-6,2.55E-6,1.43E-6,9.59E-7,5.65E-7,3.68E-7,2.74E-7,2.53E-7,2.53E-7]
MLE110 = [5.66E-4,3.31E-5,6.81E-6,2.68E-6,1.48E-6,9.84E-7,5.75E-7,3.72E-7,2.77E-7,2.54E-7,2.54E-7]
MLE120 = [9.74E-4,4.13E-5,7.55E-6,2.85E-6,1.55E-6,1.02E-6,5.87E-7,3.78E-7,2.80E-7,2.57E-7,2.56E-7]

MLE544 = [4.98E-3,1.80E-4,2.23E-5,5.91E-6,2.48E-6,1.39E-6,6.96E-7,4.03E-7,2.81E-7,2.54E-7,2.52E-7]

Masses = [5,6,7,8,9,10,12,15,20,25,30]


plt.semilogy()
plt.plot(Eric[0],Eric[1],'-og',label='EDELWEISS low mass')
plt.plot(Masses,MLE100,'-r',label=r'$100\%\,E_{thresh}$', linewidth=2)
plt.plot(Masses,MLE110,'-m',label=r'$110\%\,E_{thresh}$', linewidth=2)
plt.plot(Masses,MLE120,'-b',label=r'$120\%\,E_{thresh}$', linewidth=2)
plt.plot(Masses,MLE544,'-g',label=r'$v_{esc} = 544$', linewidth=2)
plt.legend()
plt.show()