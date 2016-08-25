import sys
import csv
from mpl_toolkits.mplot3d import axes3d
import matplotlib.cm as cm
import numpy as np
import matplotlib.pyplot as plt
import math as m

def importData(fileName):

  dataArray=[]
  file=open(fileName,'rb')
  data = csv.reader(file,delimiter=",")
  infoArray =data.next()
  for row in data:
    dataArray.append(row[1:]) 
  file.close()
  
  return infoArray,dataArray


#################

#Q = int(sys.argv[1]);
info,data=importData("/data/trafficFlowData/fft1024_track32768_runs1e1_density.05_vel9_sampleRate1.csv")
#data=importData("test.csv")
z=np.array(data).astype(float)

trackLength = int(info[5])
fftSamples=int(info[11])
numCars=int(info[1])
maxVel=int(info[3])
infoText="Track Length: {}, ".format(trackLength)
infoText=infoText+"Num of samples in each FFT: {}\n".format(fftSamples)
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
#Create S(q) plot (sum over all w for each q)
#########################################################

if 0:
  Z = np.sum(z,axis=0)
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
  qVals = range(0,3200,2)
  for q in qVals:

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
      continue

    #Save to array of full widths
    FWHM[q] = fwhm

    ################################################################
    #Find FWHM using a quadratic fit of the inverse of the peak
    ################################################################
    fitDataX1 = int(lowSideIndex+.5-3)
    fitDataX2 = int(highSideIndex-.5+3)
    numPoints = fitDataX2 - fitDataX1 + 1
    
    if (fitDataX1 < 0) or (fitDataX2 >= len(zq)):
      continue

    fitData=np.empty(zqLen)
    fitData.fill(np.NAN)
    for i in xrange(fitDataX1,fitDataX2+1):
      fitData[i] = zq[i]
    fitRange = range(fitDataX1,fitDataX2+1)
    fitParam = np.polyfit(fitRange,1.0/fitData[fitDataX1:fitDataX2+1],2)
    a = fitParam[0]
    b = fitParam[1]
    c = fitParam[2]
    
    fitDensity = .2
    fitLineX=np.empty((numPoints+4)/fitDensity+1)
    fitLineY=np.empty((numPoints+4)/fitDensity+1)
    fitLineX.fill(np.NAN)
    fitLineY.fill(np.NAN)

    j=0
    for i in np.arange(fitDataX1-2,fitDataX2+2,fitDensity):
      fitLineX[j] = i
      fitLineY[j] = a*i*i + b*i + c
      j += 1

    fitHalfMax = np.nanmax(1.0/fitLineY)/2.0

    #root = b**2.0 - 4.0*a*(c-fitHalfMax)
    root = -b**2.0 + 4.0*a*c
    if root < 0:
      print "Imaginary answer for q={}".format(q)
      print "a,b,c = ",a,b,c
      print "Max: ",np.nanmax(fitLineY)
      print "root: ",root
      continue

    #halfWidthX1 = (-b*1.0 + m.sqrt(root))/(2.0*a)
    #halfWidthX2 = (-b*1.0 - m.sqrt(root))/(2.0*a)
    #fitfwhm = halfWidthX2 - halfWidthX1
    fitPeakX = (-b)/(2.0*a)
    fitfwhm = m.sqrt(root)/a
    fitFWHM[q] = fitfwhm
#    halfWidthX1 = fitPeakX - fitfwhm/2.0
#    halfWidthX2 = fitPeakX + fitfwhm/2.0
#    print fitPeakX,fitfwhm
#
#
#    xmin=fitDataX1 - 2*fwhm
#    xmax=fitDataX2 + 2*fwhm
#    ymin=0
#    ymax=1.1*maxVal
#    fig = plt.figure()
#    plt.plot(zq);
#    plt.plot(fitData,'ro')
#    plt.plot(fitLineX,1.0/fitLineY,'g*')
#    plt.plot([lowSideIndex,highSideIndex],[halfMax,halfMax],'r-',linewidth=2.0)
#    plt.plot([halfWidthX1,halfWidthX2],[fitHalfMax,fitHalfMax],'g-',linewidth=2.0)
#
#    plt.title( 'Omega distribution for q = {}. Area={}'.format(q,np.sum(zq)) )
#    plt.xlabel('Omega')
#    plt.ylabel("Intensity")
#
#    plt.axis([xmin, xmax, ymin, ymax])
#    plt.grid(True)
#    #plt.tight_layout()
#
#    plt.savefig("/data/plots3/Q{0:04d}.png".format(q))
#    #plt.show()
#    plt.close()
  
  #####################################################################
  # Fit line to FWHM values
  #####################################################################

  xValsToFit = qVals[200:1600]
  yValsToFit = np.empty(len(xValsToFit))
  j=0
  for i in xValsToFit:
    yValsToFit[j] = fitFWHM[i]
    j+=1


  linFitParam = np.polyfit(xValsToFit,yValsToFit,1)
  m = linFitParam[0]
  b = linFitParam[1]
  print m,b

  linFitYvals=np.empty(len(xValsToFit))
  linFitYvals.fill(np.NAN)
  
  j=0
  for i in xValsToFit:
    linFitYvals[j] = m*i + b
    j+=1

  fig = plt.figure()
  plt.plot(FWHM,'r.')
  plt.plot(fitFWHM,'g.')
  plt.plot(xValsToFit,linFitYvals,'g-')

  plt.title('FWHM for the peak in each q slice\n'+infoText)
  plt.xlabel('Wave vector, q')
  plt.ylabel("FWHM")

  plt.grid(True)
  ymax = np.nanmax(fitFWHM)
  plt.axis([np.amin(qVals),np.amax(qVals),0,1.1*ymax])
  plt.tight_layout()
  plt.show()
  plt.close()


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





