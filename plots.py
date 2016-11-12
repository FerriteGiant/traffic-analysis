import sys
import csv
from mpl_toolkits.mplot3d import axes3d
import matplotlib.cm as cm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import math as m
from scipy.optimize import least_squares

def importData(fileName):

  dataArray=[]
  file=open(fileName,'rb')
  data = csv.reader(file,delimiter=",")
  infoArray =data.next()
  for row in data:
    dataArray.append(row[1:]) #[1:] skips the first entry on each row
  file.close()
  
  return infoArray,dataArray

#Pass in a parameter array 
#p[0] = height of the peak (ish...)
#p[1] = center of peak
#p[2] = FWHM
def lorentz(p,x):
  return p[0]*.5*p[2]/((x-p[1])**2 + (.5*p[2])**2)

def residuals(p,x,y):
  return lorentz(p,x) - y

def convertW(Wn,samples):
  return 2*np.pi*(Wn)/(samples*2.0)

#################

#Q = int(sys.argv[1]);
dataFiles = [
#"/data/trafficFlowData/vel5_density0.01_fft2048_track8192_runs3.01e+05_numCars81_numDataPoints5.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.03_fft2048_track8192_runs9.96e+04_numCars245_numDataPoints5.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.04_fft2048_track8192_runs1.49e+03_numCars327_numDataPoints1.0e+09_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.05_fft2048_track8192_runs5.97e+04_numCars409_numDataPoints5.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.06_fft2048_track8192_runs9.94e+02_numCars491_numDataPoints1.0e+09_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.07_fft2048_track8192_runs4.26e+04_numCars573_numDataPoints5.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.08_fft2048_track8192_runs7.45e+02_numCars655_numDataPoints1.0e+09_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.09_fft2048_track8192_runs3.31e+04_numCars737_numDataPoints5.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.10_fft2048_track8192_runs2.98e+03_numCars819_numDataPoints5.0e+09_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.11_fft2048_track8192_runs2.71e+04_numCars901_numDataPoints5.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.12_fft2048_track8192_runs4.96e+02_numCars983_numDataPoints1.0e+09_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.13_fft2048_track8192_runs2.29e+04_numCars1064_numDataPoints5.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.135_fft2048_track8192_runs4.42e+03_numCars1105_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.138_fft2048_track8192_runs4.32e+03_numCars1130_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.14_fft2048_track8192_runs8.52e+02_numCars1146_numDataPoints2.0e+09_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.145_fft2048_track8192_runs4.11e+03_numCars1187_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.15_fft2048_track8192_runs7.95e+02_numCars1228_numDataPoints2.0e+09_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.20_fft2048_track8192_runs5.96e+02_numCars1638_numDataPoints2.0e+09_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.30_fft2048_track8192_runs9.93e+02_numCars2457_numDataPoints5.0e+09_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.40_fft2048_track8192_runs7.45e+02_numCars3276_numDataPoints5.0e+09_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.50_fft2048_track8192_runs5.96e+02_numCars4096_numDataPoints5.0e+09_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.60_fft2048_track8192_runs4.96e+02_numCars4915_numDataPoints5.0e+09_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.70_fft2048_track8192_runs4.25e+02_numCars5734_numDataPoints5.0e+09_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.80_fft2048_track8192_runs3.72e+02_numCars6553_numDataPoints5.0e+09_halfdata_Hamm.csv",
#"/data/trafficFlowData/old/vel5_density0.05_fft2048_track8192_runs1e+04_numCars409_halfdata_Welch.csv",
#####################
####### Vel 9 #######
#####################
#"/data/trafficFlowData/vel9_density0.01_fft2048_track8192_runs3.01e+05_numCars81_numDataPoints5.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel9_density0.02_fft2048_track8192_runs1.50e+05_numCars163_numDataPoints5.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel9_density0.03_fft2048_track8192_runs9.96e+04_numCars245_numDataPoints5.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel9_density0.04_fft2048_track8192_runs7.47e+04_numCars327_numDataPoints5.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel9_density0.05_fft2048_track8192_runs5.97e+04_numCars409_numDataPoints5.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel9_density0.06_fft2048_track8192_runs4.97e+04_numCars491_numDataPoints5.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel9_density0.07_fft2048_track8192_runs4.26e+04_numCars573_numDataPoints5.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel9_density0.075_fft2048_track8192_runs7.95e+03_numCars614_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel9_density0.08_fft2048_track8192_runs3.73e+04_numCars655_numDataPoints5.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel9_density0.082_fft2048_track8192_runs7.28e+03_numCars671_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel9_density0.083_fft2048_track8192_runs7.19e+03_numCars679_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel9_density0.0835_fft2048_track8192_runs7.14e+03_numCars684_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel9_density0.084_fft2048_track8192_runs7.10e+03_numCars688_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel9_density0.085_fft2048_track8192_runs7.02e+03_numCars696_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel9_density0.087_fft2048_track8192_runs6.86e+03_numCars712_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel9_density0.088_fft2048_track8192_runs6.78e+03_numCars720_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel9_density0.089_fft2048_track8192_runs6.70e+03_numCars729_numDataPoints1.0e+10_halfdata_Hamm.csv",
"/data/trafficFlowData/vel9_density0.09_fft2048_track8192_runs6.62e+03_numCars737_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel9_density0.095_fft2048_track8192_runs6.28e+03_numCars778_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel9_density0.10_fft2048_track8192_runs5.96e+03_numCars819_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel9_density0.15_fft2048_track8192_runs3.98e+03_numCars1228_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel9_density0.20_fft2048_track8192_runs2.98e+03_numCars1638_numDataPoints1.0e+10_halfdata_Hamm.csv",
##"/data/trafficFlowData/vel9_density0.30_fft2048_track8192_runs1.99e+03_numCars2457_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel9_density0.50_fft2048_track8192_runs1.19e+03_numCars4096_numDataPoints1.0e+10_halfdata_Hamm.csv",

]

for file in dataFiles:

  info,data=importData(file)
  #info,data=importData("/data/trafficFlowData/fft1024_track8192_runs1e4_density.08_vel9_sampleRate1_labels.csv")
  #data=importData("test.csv")
  z=np.array(data).astype(float)
  
  trackLength = int(info[5])
  fftSamples=int(info[11])
  fftRuns = int(info[13])
  numCars=int(info[1])
  maxVel=int(info[3])
  density = float(numCars)/float(trackLength)
  densityString = "Density: {:.4f}".format(density)
  numDataPoints = int(info[15])
  infoText="Max Vel: {}, ".format(maxVel)
  infoText=infoText+"Car density: {:.03f}, ".format(density)
  infoText=infoText+"Track Length: {}, ".format(trackLength)
  infoText=infoText+"Num of time steps in each FFT: {}".format(fftSamples)
  #infoText=infoText+"Num of FFTs: {}\n".format(fftRuns)
  #infoText=infoText+"Number of cars: {}, ".format(numCars)
  #infoText=infoText+"Num of datums: {:.03e}".format(numDataPoints)
  print infoText
  print densityString
  
  
  #########################################################
  #Create S(q,w) intensity plot
  #########################################################
  
  if 1:
    ax1 = plt.subplot2grid((3,3), (0,0), colspan=3)
    ax2 = plt.subplot2grid((3,3), (1,0), colspan=3, rowspan=2)
  #  fig = plt.figure()
    #ax2.matshow(np.log(z),cmap='autumn')
    lnz =  np.log(z)
    #lnzfull = np.hstack((lnz,np.flipud(np.fliplr(lnz))))
    ax2.pcolormesh(lnz,cmap='autumn',vmin=-7,vmax=np.max(lnz))
    ax2.axvline(4096*2/(maxVel+.9),color='k')
    ax2.set_xlim(0,trackLength/(maxVel-2))
    ax2.set_ylim(0,fftSamples*2)
    ax2.grid(True)
    ax2.set_xlabel("Integer Wave Vector, q")
    ax2.set_ylabel(r"Integer frequency, $\omega_n$")
    #ax2.colorbar()
    #plt.show()
  
   # #nx: time (w)
   # #ny: space (q)
   # nx,ny = z.shape
   # Z = np.transpose(z);
   # x=np.linspace(0,nx-1,nx)
   # y=np.linspace(0,ny-1,ny)
  
   # print z.shape
   # print Z.shape
   # print len(x)
   # print len(y)
  
  
   # xmin = x.min()
   # xmax = x.max()
   # ymin = y.min()
   # ymax = y.max()
   # 
   # plt.subplots_adjust(hspace=0.5)
   # plt.subplot(121)
   # plt.hexbin(x, y, C=z, cmap=plt.cm.YlOrRd_r)
   # plt.axis([xmin, xmax, ymin, ymax])
   # plt.title("Hexagon binning")
   # cb = plt.colorbar()
   # cb.set_label('counts')
   # 
   # plt.subplot(122)
   # plt.hexbin(x, y, Z, bins='log', cmap=plt.cm.YlOrRd_r)
   # plt.axis([xmin, xmax, ymin, ymax])
   # plt.title("With a log color scale")
   # cb = plt.colorbar()
   # cb.set_label('log10(N)')
   # 
   # plt.show()
  

  if 0:
    #xmin = 0
    #xmax = len(Z)
    #ymin = 0
    #ymax = np.max(Z)*1.1
  
    fig = plt.figure()
    plt.plot(z[:,10])
  
#    plt.title("Frequency Power Spectrum\n"+infoText)
    plt.xlabel('Freq, w')
    plt.ylabel("Power")
  
#   plt.axis([xmin, xmax, ymin, ymax])
    plt.grid(True)
    plt.tight_layout()
  
    plt.show()
  
  #########################################################
  #Create S(q) plot (sum over all w for each q)
  #########################################################
  
  if 1:
    Z = np.sum(z,axis=0) #sum along time freq axis
    xmin = 0
    xmax = trackLength/(maxVel-2)
    ymin = 0
    ymax = np.max(Z)*1.1
  
    Zfull = np.hstack((Z,np.flipud(Z)))
    #ax1 = plt.subplot2grid((3,3), (0,0), colspan=3)
    ax1.plot(Z)
    ax1.axvline(4096*2/(maxVel+.9),color='k')
  
    ax1.set_title("Structure Factor\n"+infoText)
    ax1.set_xlabel('Integer Wave vector, q')
    ax1.set_ylabel(r"Sum over $\omega$")
  
    ax1.set_xlim([xmin, xmax])
    ax1.set_ylim([ymin, 450])
    ax1.grid(True)
    plt.tight_layout()
  
    plt.show()
  
  
  #########################################################
  #Plot S(q,w) in 3D
  #########################################################
  
  if 0:
    #nx: time (w)
    #ny: space (q)
    nx,ny = z.shape
    x=np.linspace(0,nx-1,nx)
    y=np.linspace(0,ny-1,ny)
    Z=np.transpose(z)
    
    #Enable to rearange data to make (0,0) in the middle
    if False:
      uLeft=z[0:nx/2+1,0:ny/2+1]
      lLeft=z[nx/2+1:nx,0:ny/2+1]
      uRight=z[0:nx/2+1,ny/2+1:ny]
      lRight=z[nx/2+1:nx,ny/2+1:ny]
      
      left=np.concatenate((lLeft,uLeft))
      right=np.concatenate((lRight,uRight))
      Z=np.concatenate((right,left),axis=1)
      
      x=np.linspace(-nx/2+1,nx/2,nx)
      y=np.linspace(-ny/2+1,ny/2,ny)
  
    #Enable to normalize as S(q,w)/S(q)
    if False:
      Z=np.empty([nx,ny])
      Sq = np.sum(z,axis=0)
      for i in range(0,ny):
        for j in range(0,nx):
          Z[j,i] = z[j,i]/Sq[i]
   
  
    X,Y=np.meshgrid(x,y)
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')  #space     #time 
    #X, Y, Z = axes3d.get_test_data(0.05)
    cset = ax.contourf(X, Y, Z, cmap=cm.coolwarm)
    ax.plot_surface(X, Y, Z, rstride=80, cstride=50, alpha=0.3,linewidth=1)
  #  cset = ax.contour(X, Y, Z, zdir='z', offset=-100, cmap=cm.coolwarm)
    cset = ax.contour(X, Y, Z, zdir='x', offset=-400, cmap=cm.coolwarm)
    cset = ax.contour(X, Y, Z, zdir='y', offset=400, cmap=cm.coolwarm)
  
    #ax.set_xlabel('X')
    #ax.set_xlim(-40, 40)
    #ax.set_ylabel('Y')
    #ax.set_ylim(-40, 40)
    #ax.set_zlabel('Z')
    #ax.set_zlim(-100, 100)
  
  #  surf=ax.plot_wireframe(np.fliplr(X),Y,np.transpose(Z),rstride=80,cstride=30,
  #                        linewidth=1,cmap=cm.coolwarm)
  
    #plt.xlabel('Omega')
    #plt.ylabel("q")
    plt.grid(True)
    #plt.tight_layout()
    #plt.savefig("data_1e5datums",ext="png", close=True)
    plt.show()
    #savefig('foo.png', bbox_inches='tight')
    #plt.close()

plt.show()

#######################################################
#######################################################
#######################################################

#vmin=0,vmax=2,
#fig.colorbar(surf, shrink=0.5, aspect=5)

#plt.figure(figsize=(16.0,9.0))
#GofR, = plt.plot(data[:,0],'ro',label="G(r)")
#PofR, = plt.plot(data[:,2],'bo',label="P(r)")
#xmin,xmax = plt.xlim()
#plt.plot([xmin,xmax],[.08,.08],'k-',lw=2)
#
##my_xticks = []
##for xtic in binDivs:
##  my_xticks.append(str(int(m.floor(xtic)))+":"+str(int((xtic-m.floor(xtic))*60)).zfill(2))
##plt.xticks(binDivs, my_xticks)
#
###plt.xticks(rotation=40)


#plt.title("")
#ax.set_zlim(0,.4e7)
#plt.legend(loc="upper right",fontsize="x-large",numpoints=1)
##plt.legend([bad,new,known],('Bad','New Dead','Known Dead'),'best',numpoints=1,
##  bbox_to_anchor=(0.9,1.05),fancybox=True,shadow=True)
##
#
##
#





