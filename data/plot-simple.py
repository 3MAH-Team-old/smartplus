# -*- coding: utf-8 -*-

#This is a very simple matplotlib script to plot results from 'results_job.txt'

from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
data_path = '/Users/chemisky/Documents/Professionnel/Methodes numériques/SMART+_free/Control/'
calcul = data_path + 'results_job.txt'

x1, y1 = np.loadtxt(calcul, usecols=(7,13), unpack=True)
plt.grid(True)

#plt.xlim(0.,1000.0)
#plt.ylim(0.,500.0)

plt.plot(x1,y1, c='blue')#,marker='b-')

#plt.xlim(0.,0.04)    
#plt.ylim(-1000,1000)  

plt.show()
