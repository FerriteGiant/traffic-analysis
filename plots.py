import sys
import csv
from mpl_toolkits.mplot3d import axes3d
import matplotlib.cm as cm
import numpy as np
import matplotlib.pyplot as plt
import math as m
from scipy.optimize import least_squares

def importData(fileName):

  dataArray=[]
  file=open(fileName,'rb')
  data = csv.reader(file,delimiter=",")
  infoArray =data.next()
  for row in data:
    dataArray.append(row[1:]) 
  file.close()
  
  return infoArray,dataArray

#Pass in a parameter array 
#p[0] = height of the peak
#p[1] = center of peak
#p[2] = FWHM
def lorentz(p,x):
  return p[0]*.5*p[2]/((x-p[1])**2 + (.5*p[2])**2)

def residuals(p,x,y):
  return lorentz(p,x) - y

#################

#Q = int(sys.argv[1]);
dataFiles = [
"/data/trafficFlowData/vel9_density0.06_fft1024_track8192_runs2e+03_numCars491_halfdata_Hamm.csv",
"/data/trafficFlowData/vel9_density0.06_fft1024_track8192_runs2e+03_numCars491_halfdata_Welch.csv",
#"/data/trafficFlowData/vel5_density0.10_fft1024_track8192_runs1e+05_numCars819_halfdata.csv",
#"/data/trafficFlowData/vel5_density0.14_fft1024_track8192_runs1e+05_numCars1146_halfdata.csv",
#"/data/trafficFlowData/vel9_density0.05_fft1024_track8192_runs1e+05_numCars409_halfdata.csv",
#"/data/trafficFlowData/vel9_density0.07_fft1024_track8192_runs1e+05_numCars573_halfdata.csv",
#"/data/trafficFlowData/vel9_density0.08_fft1024_track8192_runs1e+05_numCars655_halfdata.csv"
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
  infoText="Track Length: {}, ".format(trackLength)
  infoText=infoText+"Num of samples in each FFT: {}, ".format(fftSamples)
  infoText=infoText+"Num of FFTs: {}\n".format(fftRuns)
  infoText=infoText+"Max Vel: {}, ".format(maxVel)
  infoText=infoText+"Number of cars: {}, ".format(numCars)
  infoText=infoText+"Car density: {:.03f}".format(1.0*numCars/trackLength)
  print infoText
  
  
  #########################################################
  #Create S(q,w) intensity plot
  #########################################################
  
  if 0:
  #  fig = plt.figure()
    plt.matshow(z)
    plt.grid(True)
    plt.xlabel('q')
    plt.ylabel("omega")
    plt.colorbar()
    plt.show()
  
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
  
  #########################################################
  # Power spectrum for w at a given q
  #########################################################
  
  if 0:
    #xmin = 0
    #xmax = len(Z)
    #ymin = 0
    #ymax = np.max(Z)*1.1
  
    fig = plt.figure()
    plt.plot(z[:,0])
  
    plt.title("Frequency Power Spectrum\n"+infoText)
    plt.xlabel('Freq, w')
    plt.ylabel("Power")
  
#   plt.axis([xmin, xmax, ymin, ymax])
    plt.grid(True)
    plt.tight_layout()
  
    plt.show()
  
  #########################################################
  #Create S(q) plot (sum over all w for each q)
  #########################################################
  
  if 0:
    Z = np.sum(z,axis=0) #sum along time freq axis
    xmin = 0
    xmax = len(Z)
    ymin = 0
    ymax = np.max(Z)*1.1
  
    fig = plt.figure()
    plt.plot(Z)
  
    plt.title("Structure Factor\n"+infoText)
    plt.xlabel('Wave vector, q')
    plt.ylabel("Sum over Omega")
  
    plt.axis([xmin, xmax, ymin, ymax])
    plt.grid(True)
    plt.tight_layout()
  
    plt.show()
  
  
  #########################################################
  #Plot q slices (all omega for a given q)
  #########################################################
  
  if 1:
    globalMax = .2#np.amax(z)
    zqLen = len(z[:,0]) #Number of omega values per q
    qLen = len(z[0,:])  #Track length (ish)
    FWHM=np.empty(qLen)
    FWHM.fill(np.NAN)
    fitFWHM=np.empty(qLen)
    fitFWHM.fill(np.NAN)
  
    #Iterate over all q
    qVals = range(0,400,1)
    for q in qVals:
      #print "q = {}".format(q)
  
      normFactor = np.sum(z[:,q])
      zq = z[:,q]/normFactor

      maxVal = np.max(zq)

      maxIndex = np.argmax(zq)
      halfMax = maxVal/2
  
  
      ################################################################
      #Find FWHM via 'brute force' method
      ################################################################
  
      #iterate forward looking for low side location of half max
      for i in xrange(0,zqLen):
        if zq[i] >= halfMax:
          lowSideIndex = i-.5
          break
      try:
        lowSideIndex
      except NameError:
        print "Low side of half width not found for q = {}".format(q)
        quit()
  
      #iterate backward looking for high side location of half max
      for i in xrange(zqLen-1,-1,-1):
        if zq[i] >= halfMax:
          highSideIndex = i+.5
          break
  
      #Error Checking
      try:
        highSideIndex
      except NameError:
        print "High side of half width not found for q = {}".format(q)
        quit()
  
      if highSideIndex <= lowSideIndex:
        print "High side lower than low side for q={}".format(q)
        continue
  
      if (lowSideIndex <= 0) or (highSideIndex >= qLen):
        #print "low or high side index out of range for q={}".format(q)
        continue
  
      #Calculate FWHM and more error checking
      fwhm = highSideIndex - lowSideIndex
      if fwhm > zqLen/2:
        print "FWHM out of bounds:"
        print q,fwhm,lowSideIndex,highSideIndex
        fwhm = zqLen-highSideIndex+lowSideIndex
        print fwhm
        continue
  
      #Save to array of full widths
      FWHM[q] = fwhm
      
      ################################################################
      #Find FWHM via a non-linear least sqrs fit to a lorenzian
      ################################################################
      fitDataX1 = int(lowSideIndex+.5-4)
      fitDataX2 = int(highSideIndex-.5+4)
      if (fitDataX1 < 0) or (fitDataX2 >= len(zq)):
        continue
      xData = range(fitDataX1,fitDataX2+1)

      yData = np.empty(len(xData))
      i = 0
      for j in xData:
        yData[i] = zq[j]
        i+=1;

      x0 = (highSideIndex + lowSideIndex)/2; #guess for peak location
      p0 = np.array([maxVal,x0,fwhm]) #initial parameter guesses
  
      res = least_squares(residuals, p0, args=(xData, yData))
      #print "p0: "
      #print p0
      #print "res.x: "
      #print res.x
      A = res.x[0]
      X0 = res.x[1]
      L = res.x[2]
      fitFWHM[q] = res.x[2]
      xFitPoints = np.arange(fitDataX1,fitDataX2,.2)
      yFitPoints = lorentz(res.x,xFitPoints)


  #    ################################################################
  #    #Find FWHM using a quadratic fit of the inverse of the peak
  #    ################################################################
  #    fitDataX1 = int(lowSideIndex+.5-3)
  #    fitDataX2 = int(highSideIndex-.5+3)
  #    numPoints = fitDataX2 - fitDataX1 + 1
  #    
  #    if (fitDataX1 < 0) or (fitDataX2 >= len(zq)):
  #      continue
  #
  #    fitData=np.empty(zqLen)
  #    fitData.fill(np.NAN)
  #    for i in xrange(fitDataX1,fitDataX2+1):
  #      fitData[i] = zq[i]
  #    fitRange = range(fitDataX1,fitDataX2+1)
  #    fitParam = np.polyfit(fitRange,1.0/fitData[fitDataX1:fitDataX2+1],2)
  #    a = fitParam[0]
  #    b = fitParam[1]
  #    c = fitParam[2]
  #    
  #    fitDensity = .2
  #    fitLineX=np.empty((numPoints+4)/fitDensity+1)
  #    fitLineY=np.empty((numPoints+4)/fitDensity+1)
  #    fitLineX.fill(np.NAN)
  #    fitLineY.fill(np.NAN)
  #
  #    j=0
  #    for i in np.arange(fitDataX1-2,fitDataX2+2,fitDensity):
  #      fitLineX[j] = i
  #      fitLineY[j] = a*i*i + b*i + c
  #      j += 1
  #
  #    fitHalfMax = np.nanmax(1.0/fitLineY)/2.0
  #
  #    #root = b**2.0 - 4.0*a*(c-fitHalfMax)
  #    root = -b**2.0 + 4.0*a*c
  #    if root < 0:
  #      print "Imaginary answer for q={}".format(q)
  #      print "a,b,c = ",a,b,c
  #      print "Max: ",np.nanmax(fitLineY)
  #      print "root: ",root
  #      continue
  #
  #    halfWidthX1 = fitDataX1-3#(-b*1.0 + m.sqrt(root))/(2.0*a)
  #    halfWidthX2 = fitDataX2+3#(-b*1.0 - m.sqrt(root))/(2.0*a)
  #    #fitfwhm = halfWidthX2 - halfWidthX1
  #    fitPeakX = (-b)/(2.0*a)
  #    fitfwhm = m.sqrt(root)/a
  #    fitFWHM[q] = fitfwhm
  #    halfWidthX1 = fitPeakX - fitfwhm/2.0
  #    halfWidthX2 = fitPeakX + fitfwhm/2.0
  #    print fitPeakX,fitfwhm
  #
  #
  
      ###############
      #Plot Q slices
      ###############
  #    xmin=fitDataX1 - 2*fwhm
  #    xmax=fitDataX2 + 2*fwhm
  #    fwhmFitX1 = res.x[1]-res.x[2]/2.0
  #    fwhmFitX2 = res.x[1]+res.x[2]/2.0
  #    #ymin=0
  #    #ymax=10*maxVal#np.nanmax(maxVal,np.nanmax(1.0/fitLineY))
  #    fig = plt.figure()
  #    plt.plot(zq);
  #    plt.plot(xData,yData,'ro')
  #    plt.plot(xFitPoints,yFitPoints,'g*')
  #    plt.plot([lowSideIndex,highSideIndex],[halfMax,halfMax],'r-',linewidth=2.0)
  #    plt.plot([fwhmFitX1,fwhmFitX2],[res.x[0]/2.0,res.x[0]/2.0],'g-',linewidth=2.0)
  #
  #    plt.title( 'Omega distribution for q = {}. Area={}'.format(q,np.sum(zq)) )
  #    plt.xlabel('Omega')
  #    plt.ylabel("Intensity")
  #
  #    ymin,ymax=plt.ylim()
  #    plt.axis([xmin, xmax, ymin, ymax])
  #    plt.grid(True)
  #    #plt.tight_layout()
  #
  #    #plt.savefig("/data/tempPlots/Q{0:04d}.png".format(q))
  #    plt.show()
  #    plt.close()
    
    #####################################################################
    # Fit line to FWHM values
    #####################################################################
  
    xValsToFit = np.array(qVals[0:300])
    yValsToFit = np.empty(len(xValsToFit))
    j=0
    for i in xValsToFit:
      yValsToFit[j] = fitFWHM[i]
      j+=1
  
    indx = np.isfinite(xValsToFit) & np.isfinite(yValsToFit)
    linFitParam = np.polyfit(xValsToFit[indx],yValsToFit[indx],1)
    M = linFitParam[0]
    b = linFitParam[1]
    print M,b
  
    linFitYvals=np.empty(len(xValsToFit))
    linFitYvals.fill(np.NAN)
    
    j=0
    for i in xValsToFit:
      linFitYvals[j] = M*i + b
      j+=1
  
  #############
  # Plot Things
  #############
    fig = plt.figure()
    plt.plot(FWHM,'r.')
    plt.plot(fitFWHM,'g.-')
    plt.plot(xValsToFit,linFitYvals,'k-')
#  
    plt.title('FWHM for the peak in each q slice\n'+infoText)
    plt.xlabel('Wave vector, q')
    plt.ylabel("FWHM")
#  
    plt.grid(True)
    ymax = np.nanmax(fitFWHM)
    plt.axis([np.amin(qVals),np.amax(qVals),0,5.1])
    plt.tight_layout()
    plt.show()
#    plt.close()
  
  
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





