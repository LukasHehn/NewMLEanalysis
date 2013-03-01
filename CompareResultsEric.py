import numpy as np
import matplotlib.pyplot as plt
MLE = [[7,8,10,12,15,20,25,30],[3.04e-05,2.25e-06,1.01e-06,5.94e-07,3.74e-07,2.71e-07,2.48e-07,2.47e-07]]
Eric = [[7,8,10,12,15,20,25,30],[7.90E-4,1.17E-4,1.35E-5,4.17E-6,1.55E-6,7.34E-7,5.47E-7,4.85E-7]]
plt.semilogy()
plt.plot(Eric[0],Eric[1],'-ob',label='Eric')
plt.plot(MLE[0],MLE[1],'-or',label='MLE')
plt.legend()
plt.show()

for i in range(len(MLE[0])):
  print MLE[0][i], MLE[1][i], Eric[1][i], Eric[1][i]/MLE[1][i]