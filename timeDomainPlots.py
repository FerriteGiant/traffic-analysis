import sys
import csv
from mpl_toolkits.mplot3d import axes3d
import matplotlib.cm as cm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import math as m
from scipy.optimize import least_squares
import matplotlib as mpl
label_size = 9
mpl.rcParams['xtick.labelsize'] = label_size 
mpl.rcParams['ytick.labelsize'] = label_size 

def importData(fileName):

  dataArray=[]
  file=open(fileName,'rb')
  data = csv.reader(file,delimiter=",")
  infoArray =data.next()
  for row in data:
    dataArray.append(row) #[1:] skips the first entry on each row
  file.close()
  
  return infoArray,dataArray


#################

ax1 = plt.subplot(421)
ax2 = plt.subplot(423)
ax3 = plt.subplot(425)
ax4 = plt.subplot(427)
ax5 = plt.subplot(422)
ax6 = plt.subplot(424)
ax7 = plt.subplot(426)
ax8 = plt.subplot(428)
axisArray = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]
#Q = int(sys.argv[1]);
dataFiles = [
### Shorter burn in (10000 steps)
#"/data/trafficFlowData/timeDomain/burnIn/phi0_vel9_density0.0830_track5000_steps100000_numCars414.csv",
#"/data/trafficFlowData/timeDomain/burnIn/phi0_vel9_density0.0840_track5000_steps100000_numCars420.csv",
#"/data/trafficFlowData/timeDomain/burnIn/phi0_vel9_density0.0850_track5000_steps100000_numCars425.csv",
#"/data/trafficFlowData/timeDomain/burnIn/phi0_vel9_density0.0870_track5000_steps100000_numCars435.csv",
#
#"/data/trafficFlowData/timeDomain/phi0_vel9_density0.0830_track5000_steps100000_numCars414.csv",
#"/data/trafficFlowData/timeDomain/phi0_vel9_density0.0840_track5000_steps100000_numCars420.csv",
#"/data/trafficFlowData/timeDomain/phi0_vel9_density0.0850_track5000_steps100000_numCars425.csv",
#"/data/trafficFlowData/timeDomain/phi0_vel9_density0.0870_track5000_steps100000_numCars435.csv",
##
#"/data/trafficFlowData/timeDomain/phi0_vel5_density0.1300_track5000_steps100000_numCars650.csv",
#"/data/trafficFlowData/timeDomain/phi0_vel5_density0.1340_track5000_steps100000_numCars670.csv",
#"/data/trafficFlowData/timeDomain/phi0_vel5_density0.1380_track5000_steps100000_numCars690.csv",
#"/data/trafficFlowData/timeDomain/phi0_vel5_density0.1420_track5000_steps100000_numCars710.csv",
#
#bad"/data/trafficFlowData/timeDomain/5k2/phi0_vel9_density0.0830_track5000_steps100000_numCars414.csv",
#"/data/trafficFlowData/timeDomain/5k2/phi0_vel9_density0.0840_track5000_steps100000_numCars420.csv",
#"/data/trafficFlowData/timeDomain/5k2/phi0_vel9_density0.0850_track5000_steps100000_numCars425.csv",
#"/data/trafficFlowData/timeDomain/5k2/phi0_vel9_density0.0870_track5000_steps100000_numCars435.csv",
#
#"/data/trafficFlowData/timeDomain/5k2/phi0_vel5_density0.1300_track5000_steps100000_numCars650.csv",
#"/data/trafficFlowData/timeDomain/5k2/phi0_vel5_density0.1340_track5000_steps100000_numCars670.csv",
#"/data/trafficFlowData/timeDomain/5k2/phi0_vel5_density0.1380_track5000_steps100000_numCars690.csv",
#"/data/trafficFlowData/timeDomain/5k2/phi0_vel5_density0.1420_track5000_steps100000_numCars710.csv",
#
#"/data/trafficFlowData/timeDomain/phi0_vel9_density0.0830_track8192_steps100000_numCars679.csv",
#"/data/trafficFlowData/timeDomain/phi0_vel9_density0.0840_track8192_steps100000_numCars688.csv",
#"/data/trafficFlowData/timeDomain/phi0_vel9_density0.0850_track8192_steps100000_numCars696.csv",
#"/data/trafficFlowData/timeDomain/phi0_vel9_density0.0870_track8192_steps100000_numCars712.csv",
##
#"/data/trafficFlowData/timeDomain/phi0_vel5_density0.1300_track8192_steps100000_numCars1064.csv",
#"/data/trafficFlowData/timeDomain/phi0_vel5_density0.1340_track8192_steps100000_numCars1097.csv",
#"/data/trafficFlowData/timeDomain/phi0_vel5_density0.1380_track8192_steps100000_numCars1130.csv",
#"/data/trafficFlowData/timeDomain/phi0_vel5_density0.1420_track8192_steps100000_numCars1163.csv",
# Less than or equal test
"/data/trafficFlowData/timeDomain/ltoreq/phi0_vel9_density0.0830_track5000_steps100000_numCars414.csv",
"/data/trafficFlowData/timeDomain/ltoreq/phi0_vel9_density0.0840_track5000_steps100000_numCars420.csv",
"/data/trafficFlowData/timeDomain/ltoreq/phi0_vel9_density0.0850_track5000_steps100000_numCars425.csv",
"/data/trafficFlowData/timeDomain/ltoreq/phi0_vel9_density0.0870_track5000_steps100000_numCars435.csv",
"/data/trafficFlowData/timeDomain/ltoreq/phi0_vel5_density0.1300_track5000_steps100000_numCars650.csv",
"/data/trafficFlowData/timeDomain/ltoreq/phi0_vel5_density0.1340_track5000_steps100000_numCars670.csv",
"/data/trafficFlowData/timeDomain/ltoreq/phi0_vel5_density0.1380_track5000_steps100000_numCars690.csv",
"/data/trafficFlowData/timeDomain/ltoreq/phi0_vel5_density0.1420_track5000_steps100000_numCars710.csv",
]

a = 0
for file in dataFiles:

  info,data=importData(file)
  jamCars=np.array(data[0]).astype(int)
  
  numCars=int(info[1])
  maxVel=int(info[3])
  trackLength = int(info[5])
  density = float(numCars)/float(trackLength)
  densityString = "Density: {:.3f}".format(density)
  numTimeSteps = int(info[7])
  infoText="Max Vel: {}, ".format(maxVel)
  #infoText=infoText+"Car density: {:.03f}, ".format(density)
  infoText=infoText+"Track Length: {}, ".format(trackLength)
  infoText=infoText+"Num of time steps: {}".format(numTimeSteps)
  #infoText=infoText+"Number of cars: {}, ".format(numCars)
  #infoText=infoText+"Num of datums: {:.03e}".format(numDataPoints)
  print infoText
  
  if 0:
    xmin = 0
    xmax = numTimeSteps
    ymin = 0
    ymax = 80#np.max(jamCars)*1.1
   
    
    axisArray[a].plot(jamCars,'-',label=densityString)
    axisArray[a].annotate(densityString, xy=(-12, -12), xycoords='axes points',
                size=10, ha='right', va='top',
                bbox=dict(boxstyle='round', fc='w'))
    axisArray[a].axes.get_yaxis().set_ticks([0,80])
    if a==3 or a==7:
      axisArray[a].axes.get_xaxis().set_ticks([0,numTimeSteps])
    else:
      axisArray[a].axes.get_xaxis().set_ticks([])
    if a==1:
      ax1.set_title('Vmax={}'.format(maxVel),fontsize=15)
    if a==5:
      ax5.set_title('Vmax={}'.format(maxVel),fontsize=15)
    ax4.set_xlabel('Time Steps')
    ax8.set_xlabel('Time Steps')
    if a<4:
      axisArray[a].set_ylabel("$N_0$")
   
    axisArray[a].set_xlim([xmin, xmax])
    axisArray[a].set_ylim([ymin, ymax])
   # ax1.grid(False)
    a = a+1
   
### Histogram
  if 1:
    xmin = 0
    xmax = numTimeSteps
    ymin = 0
    ymax = 80#np.max(jamCars)*1.1
   
    
    jamFrac = jamCars.astype(float)/numCars
    bins = np.linspace(0,.12,26)
    binWidth = bins[1]-bins[0]
    jamFrac = jamCars.astype(float)/numCars
    results,edges=np.histogram(jamFrac,bins,normed=1)
    axisArray[a].bar(edges[:-1],(results*binWidth),binWidth)
    #axisArray[a].plt(jamCars,'-',label=densityString)
    #axisArray[a].legend(loc="upper right",ncol=1,shadow=True,fancybox=True,fontsize="small",numpoints=None)
    axisArray[a].annotate(densityString, xy=(-12, -12), xycoords='axes points',
                size=10, ha='right', va='top',
                bbox=dict(boxstyle='round', fc='w'))
    axisArray[a].axes.get_yaxis().set_ticks([0,.25,.5])
    axisArray[a].axes.get_xaxis().set_ticks([.06,.12])
  #  if a==3 or a==7:
  #    axisArray[a].axes.get_xaxis().set_ticks([0,numTimeSteps])
  #  else:
  #    axisArray[a].axes.get_xaxis().set_ticks([])
    if a==1:
      ax1.set_title('Vmax={}'.format(maxVel),fontsize=15)
    if a==5:
      ax5.set_title('Vmax={}'.format(maxVel),fontsize=15)
    ax4.set_xlabel(r'$\phi_0/N$')
    ax8.set_xlabel(r'$\phi_0/N$')
  #  if a<4:
  #    axisArray[a].set_ylabel("$N_0$")
   
    axisArray[a].set_xlim([xmin, bins[-1]])
    axisArray[a].set_ylim([ymin, .5])
    axisArray[a].grid(True)
    a = a+1
   
plt.subplots_adjust(hspace=.25)  
#plt.tight_layout() 
#plt.suptitle("Number of Cars in a Jam",fontsize=16)
plt.suptitle("Probability distribution for finding some fraction of cars in a Jam\nTracksize = {}".format(trackLength),fontsize=16)
plt.show()
 
