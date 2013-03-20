import numpy as np
import matplotlib.pyplot as plt

EricID3 = [[7,8,10,12,15,20,25,30],[7.90E-4,1.17E-4,1.35E-5,4.17E-6,1.55E-6,7.34E-7,5.47E-7,4.85E-7]]

Xenon100 = [[5.8,6,7,10,20],[1.00E-003,5.00E-004,2.00E-005,2.00E-007,7.00E-009]]

CDMS = [[5,6,7,8,9,10,11,12],[3.00E-004,1.10E-004,6.00E-005,3.50E-005,2.60E-005,2.00E-005,1.60E-005,1.50E-005]]

Eigene1 = [[8,10,15,30],[6.67E-5, 8.55E-6, 1.14E-6, 3.97E-7]]
Eigene2 = [[8,10,15,30],[4.70E-5, 7.38E-6, 1.09E-6, 3.90E-7]]


plt.semilogy()
plt.plot(Eigene1[0],Eigene1[1],'or-',label='with own gamma-cut', linewidth=1)
plt.plot(Eigene2[0],Eigene2[1],'ob-',label='without gamma-cut', linewidth=1)
plt.plot(EricID3[0],EricID3[1],label='EDELWEISS Low Threshold (ID3)', linewidth=2,ls='-',color='darkslategray')
plt.plot(CDMS[0],CDMS[1],label='CDMS Low Threshold', linewidth=2,ls=':',color='darkslategray')
plt.plot(Xenon100[0],Xenon100[1],label='XENON100 (2012)', linewidth=2,ls='--',color='darkslategray')
plt.legend()
plt.xlabel('WIMP mass [GeV]',fontsize=18)
plt.ylabel('WIMP-nucleon cross section [pb]',fontsize=18)
plt.ylim(1e-7,1e-3)
plt.xlim(5,31)
plt.grid(True)
#plt.grid(True,which="both",ls="-",color='black')
#plt.tick_params(axis='both', which='major', labelsize=18)
#plt.tick_params(axis='both', which='minor', labelsize=18)
plt.show()