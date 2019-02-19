#example of plotting a complete DIA along a 2D surface
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.mlab as ml

#import data from autoDIA output file
a = np.genfromtxt('Scan2D_DIA.txt',skip_header=4)

#create plot and subplots and layout
f, axes = plt.subplots(2, 3,figsize=(5.25, 5),dpi=300,gridspec_kw = {'width_ratios':[5,5, 0.5]})
f.delaxes(axes[0,1])
f.delaxes(axes[0,2])
axes[0,0].set_aspect(1, adjustable='box')
axes[1,0].set_aspect(1, adjustable='box')
axes[1,1].set_aspect(1, adjustable='box')
axes[1,2].set_aspect(10, adjustable='box')
f.subplots_adjust(wspace=0, hspace=0)

"""The second column (keep in mind that Numpy numbering starts with 0) represents the first bond distance,
the third column contains information about the second bond distance, the 6th column contains the total energy
"""
x = a[:,1]
y = a[:,2]
z = a[:,5]

#getting the limits for the axis
xmin = np.amin(x)
xmax = np.amax(x)
ymin = np.amin(y)
ymax = np.amax(y)

#Manually select lower and upper energy limits
Emin=-180
Emax=150

#create array for contour lines which should be plottet. Spacing 10 kcal/mol
b = np.linspace(Emin,Emax,(Emax-Emin)//10)

#create the griddata with z values needed for contour and pcolormesh plotting
xi = np.linspace(xmin, xmax, len(x))
yi = np.linspace(ymin, ymax, len(y))
zi = ml.griddata(x, y, z, xi, yi, interp='linear')

#plotting contour lines and colormesh for the total energy subplot
axes[0,0].contour(xi, yi, zi, b, linewidths = 0.5, colors = 'k',linestyles='solid')
z1_plot = axes[0,0].pcolormesh(xi, yi, zi,vmin=Emin,vmax=Emax, cmap = plt.get_cmap('rainbow'))

#setting axes
axes[0,0].set_title('Electronic energy')
axes[0,0].set_ylabel('Bond distance 1')
axes[0,0].set_yticks([1.5,2.0,2.5,3.0])

#plotting colorbar
cb = plt.colorbar(z1_plot,cax=axes[1,2], ticks=[-150,-100, -50, 0, 50, 100, 150]) 
cb.ax.set_yticklabels([-150,-100, -50, 0, 50, 100, 150],ha='right')
cb.ax.yaxis.set_tick_params(pad=25)
axes[0,0].set_xticklabels([])

#column 7 contains information about the interaction energy
z = a[:,6]

#create the griddata with z values needed for contour and pcolormesh plotting
zi = ml.griddata(x, y, z, xi, yi, interp='linear')

#plotting contour lines and colormesh for the interaction subplot
axes[1,0].contour(xi, yi, zi, b, linewidths = 0.5, colors = 'k',linestyles='solid')
axes[1,0].pcolormesh(xi, yi, zi,vmin=Emin,vmax=Emax, cmap = plt.get_cmap('rainbow'))

#setting axes
axes[1,0].set_title('Interaction energy')
axes[1,0].set_ylabel('Bond distance 1')
axes[1,0].set_xlabel('Bond distance 2')
axes[1,0].set_yticks([1.5,2.0,2.5,3.0])

#column 8 contains information about the distortion energy
z = a[:,7]

#create the griddata with z values needed for contour and pcolormesh plotting
zi = ml.griddata(x, y, z, xi, yi, interp='linear')

#plotting contour lines and colormesh for the distortion subplot
axes[1,1].contour(xi, yi, zi, b, linewidths = 0.5, colors = 'k',linestyles='solid')
axes[1,1].pcolormesh(xi, yi, zi,vmin=Emin,vmax=Emax, cmap = plt.get_cmap('rainbow'))

#setting axes
axes[1,1].set_title('Distortion energy')
axes[1,1].set_xlabel('Bond distance 2')
axes[1,1].set_yticklabels([])

#layout and saving the figure
f.tight_layout()
plt.savefig('2D.png')
