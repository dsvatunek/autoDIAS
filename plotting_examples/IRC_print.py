#example of plotting a complete DIA along an IRC as shown in Figure 6
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

#create plot
f= plt.figure(figsize=(2.5, 2.5),dpi=300)
ax1 = f.add_subplot(1,1,1)

#read data from autoDIA output
a = np.genfromtxt('IRC1_DIA.txt',skip_header=4)

#set axis
ax1.set_ylabel('Energy')
ax1.set_xlabel('Forming bond distance')
ax1.axis([np.amax(a[:,2])+0.1, np.amin(a[:,2])-0.1, np.amin(a[:,6])-20, np.amax(a[:,7])+20])
ax1.axhline(y=0, color='grey', linestyle='--',linewidth = 0.5)

#the second geometric parameter (column 3 - keep in mind that Numpy counting starts at 0) representing a bond length is used as X-Value
x = a[:,2]
#In this case column 7 contains information about interaction energy
y = a[:,6]
ax1.plot(x,y,'-',color='g',linewidth = 1)

#In this case column 8 contains information about total distortion energy
y = a[:,7]
ax1.plot(x,y,'-',color='r',linewidth = 1)

#In this case column 9 contains information about the acrolein distortion energy
y = a[:,8]
ax1.plot(x,y,'-',color='lightcoral',linewidth = 0.5)

#In this case column 10 contains information about the butadiene distortion energy
y = a[:,9]
ax1.plot(x,y,'-',color='orange',linewidth = 0.5)

#In this case column 6 contains information about the total energy. Plottet last to be on top.
y = a[:,5]
ax1.plot(x,y,'-',color='k',linewidth = 1)

#creating the legend
black_patch = mpatches.Patch(color='k', label='Total energy')
red_patch = mpatches.Patch(color='r', label='Distortion energy')
green_patch = mpatches.Patch(color='g', label='Interaction energy')
coral_patch = mpatches.Patch(color='lightcoral', label='Acrolein distortion energy')
orange_patch = mpatches.Patch(color='orange', label='Butadiene distortion energy')
plt.legend(handles=[black_patch,red_patch,coral_patch,orange_patch, green_patch],prop={'size': 5},loc=3,frameon=False)

f.tight_layout()

f.savefig('IRC.png')
