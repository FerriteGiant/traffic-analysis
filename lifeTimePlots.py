import sys
import csv
from mpl_toolkits.mplot3d import axes3d
import matplotlib.cm as cm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import math as m
import itertools
from scipy.optimize import least_squares
import peakutils
import scipy.signal as sig
import matplotlib as mpl
label_size = 14
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size

def importData(fileName):

  dataArray=[]
  file=open(fileName,'rb')
  data = csv.reader(file,delimiter=",")
  infoArray =data.next()
  maxVel=int(infoArray[3])
  trackLength = int(infoArray[5])
  maxQ = int(trackLength/maxVel)
  for row in data:
    dataArray.append(row[1:maxQ]) #[1:] skips the first entry on each row
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

def convertQ(Wn,trackLength):
  return 2*np.pi*(Wn)/(trackLength)

def two_largest(numbers):
    x = 0
    thresholdFrac = .025
    m1 = m2 = float('-inf')
    x1 = x2 = -1
    for m in numbers:
        if m > m2:
            if m >= m1:
                m1, m2 = m, m1            
                x1, x2 = x, x1
            else:
                m2 = m
                x2 = x
        x = x + 1
    if m1 < m2*thresholdFrac:
      return [x2]
    if m2 < m1*thresholdFrac:
      return [x1]
    if x1<x2:
      return [x1,x2]
    else:
      return [x2,x1]

markers = itertools.cycle(('o', 'v', '8', 's', '^','p', '*', 'h', '<', 'H', 'D', '>','d'))
colors = itertools.cycle(('b', 'g', 'r', 'brown','c', 'm', 'y', 'k','darksalmon'))
#markers.next()
#colors.next()
#################

HWplt = plt.figure(figsize=(10,6))
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
#"/data/trafficFlowData/vel5_density0.140_fft2048_track8192_runs2.13e+04_numCars1146_numDataPoints5.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.145_fft2048_track8192_runs2.06e+04_numCars1187_numDataPoints5.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel5_density0.150_fft2048_track8192_runs1.99e+04_numCars1228_numDataPoints5.0e+10_halfdata_Hamm.csv",
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
"/data/trafficFlowData/vel9_density0.08_fft2048_track8192_runs3.73e+04_numCars655_numDataPoints5.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel9_density0.082_fft2048_track8192_runs7.28e+03_numCars671_numDataPoints1.0e+10_halfdata_Hamm.csv",
"/data/trafficFlowData/vel9_density0.083_fft2048_track8192_runs7.19e+03_numCars679_numDataPoints1.0e+10_halfdata_Hamm.csv",
"/data/trafficFlowData/vel9_density0.0835_fft2048_track8192_runs7.14e+03_numCars684_numDataPoints1.0e+10_halfdata_Hamm.csv",
"/data/trafficFlowData/vel9_density0.084_fft2048_track8192_runs7.10e+03_numCars688_numDataPoints1.0e+10_halfdata_Hamm.csv",
"/data/trafficFlowData/vel9_density0.085_fft2048_track8192_runs7.02e+03_numCars696_numDataPoints1.0e+10_halfdata_Hamm.csv",
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
#####################
####### Vel 12 ######
#####################
#"/data/trafficFlowData/vel12_density0.020_fft2048_track8192_runs8.99e+03_numCars163_numDataPoints3.0e+09_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel12_density0.030_fft2048_track8192_runs1.99e+04_numCars245_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel12_density0.040_fft2048_track8192_runs1.49e+04_numCars327_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel12_density0.045_fft2048_track8192_runs1.06e+04_numCars368_numDataPoints8.0e+09_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel12_density0.050_fft2048_track8192_runs1.19e+04_numCars409_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel12_density0.060_fft2048_track8192_runs9.94e+03_numCars491_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel12_density0.070_fft2048_track8192_runs8.52e+03_numCars573_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel12_density0.075_fft2048_track8192_runs7.95e+03_numCars614_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel12_density0.080_fft2048_track8192_runs7.45e+03_numCars655_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel12_density0.090_fft2048_track8192_runs6.62e+03_numCars737_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel12_density0.100_fft2048_track8192_runs5.96e+03_numCars819_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel12_density0.450_fft2048_track8192_runs1.06e+03_numCars3686_numDataPoints8.0e+09_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel12_density0.550_fft2048_track8192_runs1.08e+03_numCars4505_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel12_density0.650_fft2048_track8192_runs9.17e+02_numCars5324_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/vel12_density0.750_fft2048_track8192_runs7.94e+02_numCars6144_numDataPoints1.0e+10_halfdata_Hamm.csv",
#####################
##### Track 4096 ####
#####################
#"/data/trafficFlowData/track4096/vel5_density0.070_fft2048_track4096_runs1.71e+04_numCars286_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/track4096/vel5_density0.130_fft2048_track4096_runs9.18e+03_numCars532_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/track4096/vel5_density0.138_fft2048_track4096_runs8.64e+03_numCars565_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/track4096/vel5_density0.145_fft2048_track4096_runs8.23e+03_numCars593_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/track4096/vel5_density0.160_fft2048_track4096_runs7.45e+03_numCars655_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/track4096/vel5_density0.200_fft2048_track4096_runs5.96e+03_numCars819_numDataPoints1.0e+10_halfdata_Hamm.csv",
##
#"/data/trafficFlowData/track4096/vel9_density0.040_fft2048_track4096_runs1.50e+04_numCars163_numDataPoints5.0e+09_halfdata_Hamm.csv",
#"/data/trafficFlowData/track4096/vel9_density0.060_fft2048_track4096_runs9.96e+03_numCars245_numDataPoints5.0e+09_halfdata_Hamm.csv",
#"/data/trafficFlowData/track4096/vel9_density0.080_fft2048_track4096_runs1.19e+04_numCars327_numDataPoints8.0e+09_halfdata_Hamm.csv",
#"/data/trafficFlowData/track4096/vel9_density0.090_fft2048_track4096_runs1.06e+04_numCars368_numDataPoints8.0e+09_halfdata_Hamm.csv",
#"/data/trafficFlowData/track4096/vel9_density0.100_fft2048_track4096_runs9.55e+03_numCars409_numDataPoints8.0e+09_halfdata_Hamm.csv",
#"/data/trafficFlowData/track4096/vel9_density0.150_fft2048_track4096_runs6.36e+03_numCars614_numDataPoints8.0e+09_halfdata_Hamm.csv",
######################
##### Track 16384 ####
######################
#"/data/trafficFlowData/track16384/vel5_density0.130_fft2048_track16384_runs2.29e+03_numCars2129_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/track16384/vel5_density0.138_fft2048_track16384_runs2.16e+03_numCars2260_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/track16384/vel5_density0.145_fft2048_track16384_runs2.06e+03_numCars2375_numDataPoints1.0e+10_halfdata_Hamm.csv",
##
#"/data/trafficFlowData/track16384/vel9_density0.060_fft2048_track16384_runs4.97e+03_numCars983_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/track16384/vel9_density0.070_fft2048_track16384_runs4.26e+03_numCars1146_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/track16384/vel9_density0.075_fft2048_track16384_runs3.98e+03_numCars1228_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/track16384/vel9_density0.080_fft2048_track16384_runs3.73e+03_numCars1310_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/track16384/vel9_density0.083_fft2048_track16384_runs3.59e+03_numCars1359_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/track16384/vel9_density0.085_fft2048_track16384_runs3.51e+03_numCars1392_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/track16384/vel9_density0.090_fft2048_track16384_runs3.31e+03_numCars1474_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/track16384/vel9_density0.095_fft2048_track16384_runs3.14e+03_numCars1556_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/track16384/vel9_density0.100_fft2048_track16384_runs2.98e+03_numCars1638_numDataPoints1.0e+10_halfdata_Hamm.csv",
################################
#### VARYING TRACK LENGTH ######
###############################
#"/data/trafficFlowData/varTrack/vel9_density0.080_fft2048_track2048_runs3.00e+04_numCars163_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.080_fft2048_track4096_runs1.49e+04_numCars327_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.080_fft2048_track8192_runs7.45e+03_numCars655_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.080_fft2048_track16384_runs3.73e+03_numCars1310_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.080_fft2048_track32768_runs1.86e+03_numCars2621_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.083_fft2048_track2048_runs2.89e+04_numCars169_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.083_fft2048_track4096_runs1.44e+04_numCars339_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.083_fft2048_track8192_runs7.19e+03_numCars679_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.083_fft2048_track16384_runs3.59e+03_numCars1359_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.083_fft2048_track32786_runs1.79e+03_numCars2721_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.084_fft2048_track2048_runs2.84e+04_numCars172_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.084_fft2048_track4096_runs1.42e+04_numCars344_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.084_fft2048_track8192_runs7.10e+03_numCars688_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.084_fft2048_track16384_runs3.55e+03_numCars1376_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.084_fft2048_track32768_runs1.77e+03_numCars2752_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.085_fft2048_track2048_runs2.81e+04_numCars174_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.085_fft2048_track4096_runs1.40e+04_numCars348_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.085_fft2048_track8192_runs7.02e+03_numCars696_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.085_fft2048_track16384_runs3.51e+03_numCars1392_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.085_fft2048_track32768_runs1.75e+03_numCars2785_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.087_fft2048_track2048_runs2.74e+04_numCars178_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.087_fft2048_track4096_runs1.37e+04_numCars356_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.087_fft2048_track8192_runs6.86e+03_numCars712_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.087_fft2048_track16384_runs3.43e+03_numCars1425_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel9_density0.087_fft2048_track32768_runs1.71e+03_numCars2850_numDataPoints1.0e+10_halfdata_Hamm.csv",
#
#"/data/trafficFlowData/varTrack/vel5_density0.130_fft2048_track2048_runs1.84e+04_numCars266_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.130_fft2048_track4096_runs9.18e+03_numCars532_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.130_fft2048_track8192_runs4.59e+03_numCars1064_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.130_fft2048_track16384_runs2.29e+03_numCars2129_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.130_fft2048_track32768_runs1.15e+03_numCars4259_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.134_fft2048_track2048_runs1.78e+04_numCars274_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.134_fft2048_track4096_runs8.91e+03_numCars548_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.134_fft2048_track8192_runs4.45e+03_numCars1097_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.134_fft2048_track16384_runs2.22e+03_numCars2195_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.134_fft2048_track32786_runs1.11e+03_numCars4393_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.138_fft2048_track2048_runs1.73e+04_numCars282_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.138_fft2048_track4096_runs8.64e+03_numCars565_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.138_fft2048_track8192_runs4.32e+03_numCars1130_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.138_fft2048_track16384_runs2.16e+03_numCars2260_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.138_fft2048_track32768_runs1.08e+03_numCars4521_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.142_fft2048_track2048_runs1.68e+04_numCars290_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.142_fft2048_track4096_runs8.40e+03_numCars581_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.142_fft2048_track8192_runs4.20e+03_numCars1163_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.142_fft2048_track16384_runs2.10e+03_numCars2326_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.142_fft2048_track32768_runs1.05e+03_numCars4653_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.146_fft2048_track2048_runs1.63e+04_numCars299_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.146_fft2048_track4096_runs8.16e+03_numCars598_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.146_fft2048_track8192_runs4.08e+03_numCars1196_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.146_fft2048_track16384_runs2.04e+03_numCars2392_numDataPoints1.0e+10_halfdata_Hamm.csv",
#"/data/trafficFlowData/varTrack/vel5_density0.146_fft2048_track32768_runs1.02e+03_numCars4784_numDataPoints1.0e+10_halfdata_Hamm.csv",
]

z=0
#dataLabelArray = []
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
  if numCars == 684:
    densityString = "Density: {:.4f}".format(density)
  else:
    densityString = "Density: {:.3f}".format(density)
  trackString = "Track: {}".format(trackLength)
  numDataPoints = int(info[15])
  infoText="Vmax: {}, ".format(maxVel)
  infoText=infoText+"Track Length: {}".format(trackLength)
  #infoText=infoText+"Num of time steps in each FFT: {}, ".format(fftSamples)
  #infoText=infoText+"Num of FFTs: {}\n".format(fftRuns)
  #infoText=infoText+"Number of cars: {}, ".format(numCars)
  #infoText=infoText+"Car density: {:.03f}".format(density)
  #infoText=infoText+"Num of datums: {:.03e}".format(numDataPoints)
  print "\n"
  print infoText
  print densityString, numCars
  print "Track: {}\n".format(trackLength)
  
  
  
  #########################################################
  #Plot q slices (all omega for a given q)
  #########################################################
  
  zqLen = len(z[:,0]) #Number of omega values per q
  qLen = len(z[0,:])  #Track length (ish)
  jamBool = False
  FWHM=np.empty(qLen)
  FWHM.fill(np.NAN)
  fitFWHM=np.empty(qLen)
  fitFWHM.fill(np.NAN)
  peakArea=np.empty(qLen)
  peakArea.fill(np.NAN)
  
  #Iterate over q
  maxQ = int(trackLength/(maxVel+1.5))
  qVals = range(10,maxQ,maxQ/50)
  for q in qVals:
    #print "q = {}".format(q)
  
    normFactor = np.sum(z[:,q])
    zq = z[:,q]/normFactor

    smoothedData = sig.savgol_filter(zq,65,1,mode='wrap')
    peakIndexes = peakutils.peak.indexes(zq, thres=0.009, min_dist=450)
    #If too many peaks are found, get the biggest two
    if peakIndexes.size > 2:
      peakValues=np.zeros(peakIndexes.size)
      p=0
      for pI in peakIndexes:
        peakValues[p] = zq[pI]
        p = p+1
      maxPeaks = two_largest(peakValues)#returns locations of max peaks in peakIndexes array 
      peakIndexes = np.array([peakIndexes[i] for i in maxPeaks])

    if peakIndexes.size > 1:
      jamBool = True #We have now seen two peaks at least once

#    if peakIndexes.size > 2:
#      if peakIndexes[2] > fftSamples*2-200: #ignore the wrap around
#        peakIndexes = peakIndexes[0:2]
#      else:
#        print("Error: Found more than two peaks. Discarding q={}.".format(q))

    if peakIndexes.size == 0:
      print("Error: No peaks found. Discarding q={}".format(q))
      continue

#XXX    if jamBool and peakIndexes.size == 1:# and peakIndexes[0]>2500:
#XXX      print("Wrong peak? Skipping q={}".format(q))
#XXX      continue
#XXX
#XXX    for peakI in [peakIndexes[0]]:
    for peakI in [peakIndexes[-1]]:

      maxVal = zq[peakI]
  
      halfMax = maxVal/2.0
      tenthMax = maxVal/10.0
      #areaFrac = maxVal/500.0
    
    
      ################################################################
      #Find FWHM via 'brute force' method
      ################################################################
      fitDataX1 = -1
      fitDataX2 = -1
      highSideIndex = -1
      lowSideIndex = -1
    
      #iterate forward from the peak looking for 
      #the high side index of half max
      for i in xrange(peakI,zqLen):
        if zq[i] <= halfMax:
          highSideIndex = i-.5
          break
      if highSideIndex == -1:
        print "High side of half width not found. Discarding q = {}".format(q)
        continue

      for i in xrange(int(highSideIndex+.5),zqLen-5):
        slope = (smoothedData[i+5]-smoothedData[i-1])/6
        fitDataX2 = i
        if zq[i] <= tenthMax or slope>0:
          break

      areaMaxX = fitDataX2 #In case we are super close to the end
      for i in xrange(int(fitDataX2),zqLen-8):
        if fitDataX2>(zqLen-9):
          areaMaxX = zqLen-1
          break
        slope = (smoothedData[i+8]-smoothedData[i-1])/9
        areaMaxX = i
        #if zq[i] < areaFrac or slope>0:
        if slope>0:
          areaMaxX = i+4
          if areaMaxX>zqLen:
            areaMaxX = zqLen
          break

      #iterate backward from the peak looking for 
      #the low side index of half max
      for i in xrange(peakI,-1,-1):
        if zq[i] <= halfMax:
          lowSideIndex = i+.5
          break
      #Error Checking
      if lowSideIndex == -1:
        print "Low side of half width not found. Discarding q = {}".format(q)
        continue

      for i in xrange(int(lowSideIndex-.5),-1,-1):
        fitDataX1 = i
        if zq[i] <= tenthMax:
          break
      
      areaMinX = fitDataX1 #In case we are super close to zero
      for i in xrange(int(fitDataX1),8,-1):
        if fitDataX1<9:
          areaMinX = 0
          break
        slope = (smoothedData[i+1]-smoothedData[i-8])/9
        areaMinX = i
        #if zq[i] <= areaFrac or slope<0:
        if slope<0:
          areaMinX = i-4
          if areaMinX<0:
            areaMinX = 0
          break
    
      if highSideIndex <= lowSideIndex:
        print "High side lower than low side for q={}".format(q)
        continue
    
      if (lowSideIndex <= 0) or (highSideIndex >= zqLen):
        print "low or high side index out of range for q={}".format(q)
        continue
      
      #Calculate FWHM and more error checking
      fwhm = highSideIndex - lowSideIndex
      if fwhm > zqLen/4:
        print "FWHM out of bounds:"
        print q,fwhm,lowSideIndex,highSideIndex
        fwhm = zqLen-highSideIndex+lowSideIndex
        print fwhm
        continue
      if peakI<fwhm:
        print "Peak too close to zero. Discarding q={}".format(q)
        continue
    
      #Save to array of full widths
      FWHM[q] = fwhm
      peakArea[q] = np.sum(zq[areaMinX:areaMaxX])

      
      ################################################################
      #Find FWHM via a non-linear least sqrs fit to a lorenzian
      ################################################################
      #fitDataX1 = int(lowSideIndex+.5-4)
      #fitDataX2 = int(highSideIndex-.5+4)
      if (fitDataX1 < 0) or (fitDataX2 >= len(zq)):
        print "fitDataX1 or X2 out of range, discarding q={}".format(q)
        continue
      xData = np.arange(fitDataX1,fitDataX2+1)
  
      yData = np.empty(len(xData))
      i = 0
      for j in xData:
        yData[i] = zq[j]
        i+=1;
  
      x0 = peakI; #guess for peak location
      p0 = np.array([maxVal,x0,fwhm]) #initial parameter guesses
    
    #    plt.plot(xData,yData,'o-')
    #    plt.show()
      res = least_squares(residuals, p0, args=(xData, yData))
      #print "p0: "
      #print p0
      #print "res.x: "
      #print res.x
      A = res.x[0]
      X0 = res.x[1]
      L = res.x[2]
      if L>fwhm*2.0:
        print "Fit FWHM much larger than guess. Discarding q={}".format(q)
        print fwhm,L
        continue
      if L<fwhm/2.0:
        print "Fit FWHM much smaller than guess. Discarding q={}".format(q)
        print fwhm,L
        continue
      #print X0
      #print A
      fitFWHM[q] = L
      xFitPoints = np.arange(fitDataX1,fitDataX2,.2)
      yFitPoints = lorentz(res.x,xFitPoints)


      ###############
      #Plot Q slices
      ###############
      if 0:
        #xmin=convertW(fitDataX1 - 4*fwhm,fftSamples)
        #xmax=convertW(fitDataX2 + 4*fwhm,fftSamples)
        xmin=convertW(areaMinX,fftSamples)
        xmax=convertW(areaMaxX,fftSamples)
        #print xmin,xmax
        fwhmFitX1 = convertW(res.x[1]-res.x[2]/2.0,fftSamples)
        fwhmFitX2 = convertW(res.x[1]+res.x[2]/2.0,fftSamples)
        #ymin=0
        #ymax=10*maxVal#np.nanmax(maxVal,np.nanmax(1.0/fitLineY))
        fig = plt.figure()
        zqX = convertW(np.arange(0,zqLen,1),fftSamples)
        plt.plot(zqX,zq,'-')
        plt.semilogy(zqX,smoothedData,'-')
        plt.semilogy(convertW(xData,fftSamples),yData,'ro')
        plt.semilogy(convertW(xFitPoints,fftSamples),yFitPoints,'g*')
        plt.semilogy([convertW(lowSideIndex,fftSamples),convertW(highSideIndex,fftSamples)],[halfMax,halfMax],'r-',linewidth=2.0)
        plt.semilogy([fwhmFitX1,fwhmFitX2],[max(yFitPoints)/2.0,max(yFitPoints)/2.0],'g-',linewidth=2.0)
        plt.axvline(xmin,color='c')
        plt.axvline(xmax,color='c')
        ymin,ymax=plt.ylim()
        for x in convertW(peakIndexes,fftSamples):
          plt.axvline(x,color='k')
        plt.title( 'Omega distribution for q = {} and density = {:.03f}. Area={:.03f}'.format(q,density,peakArea[q]) )
        plt.xlabel('Omega')
        plt.ylabel("Intensity")
  
        plt.axis([0, 6.3, 0, ymax])
        #plt.axis([xmin, xmax, 0, ymax])
        plt.grid(True)
        plt.tight_layout()
  
        plt.savefig("/data/Qslices/peak_V{:d}_Q{:04d}_d{:.03f}_2048_hamm.png".format(maxVel,q,density))
        #plt.show()
        plt.close()
  
    #####################################################################
    # Fit line to FWHM values
    #####################################################################
  
   # xValsToFit = np.array(qVals[0:80])
   # yValsToFit = np.empty(len(xValsToFit))
   # j=0
   # for i in xValsToFit:
   #   yValsToFit[j] = fitFWHM[i]
   #   j+=1
  
   # indx = np.isfinite(xValsToFit) & np.isfinite(yValsToFit)
   # linFitParam = np.polyfit(xValsToFit[indx],yValsToFit[indx],1)
   # M = linFitParam[0]
   # b = linFitParam[1]
   # print M,b
  
   # linFitYvals=np.empty(len(xValsToFit))
   # linFitYvals.fill(np.NAN)
   # 
   # j=0
   # for i in xValsToFit:
   #   linFitYvals[j] = M*i + b
   #   j+=1
  
  #############
  # Plot Things
  #############
  qDataPI = np.array( [convertQ(qq,trackLength) for qq in range(0,qLen)] )
  maxQrad = convertQ(maxQ,trackLength)
  ax1 = plt.subplot(211)
  m = markers.next()
  c = colors.next()
  ax1.plot(qDataPI,convertW(fitFWHM/2.0,fftSamples),color=c,linestyle='',marker=m,label=densityString)
  #ax1.plot(qDataPI,convertW(fitFWHM/2.0,fftSamples),color=c,linestyle='',marker=m,label=trackString)
 # plt.plot(gamma,convertW(fitFWHM,fftSamples)/2.0,'g.')
 # plt.plot(gamma,gamma,'b')
  #plt.plot(xValsToFit,linFitYvals,'k-')

  ax1.set_title('Free peak in each q slice\n'+infoText)

  ax1.set_ylabel('HWHM (rad / time step)',fontsize=16)
  #ymax = np.nanmax(convertW(fitFWHM,fftSamples))/2*1.1
  ax1.grid(b=True, which='major', color='k', linestyle='-')
  #plt.grid(b=True, which='minor', color='r', linestyle='-', alpha=0.2)
  ax1.set_xlim([0,maxQrad])
  ax1.set_ylim([0,.08])
  ax1.set_xticklabels([])

  ax2 = plt.subplot(212)
  ax2.plot(qDataPI,peakArea,'k',color=c,linestyle='',marker=m,label=densityString)
  ax2.legend(loc="lower right",ncol=2,shadow=True,fancybox=True,fontsize=13,numpoints=1)
  ax2.set_xlim([0,maxQrad])
  ax2.set_ylim([0,1.1])
  ax2.set_ylabel("Peak area",fontsize=16)
  #ticks = ax2.xaxis.get_majorticklocs()
  #tickLabels = ["$\pi$/{:.0f}".format(np.pi/i) for i in ticks]
  #tickLabels[0] = "0"
  #ax2.set_xticklabels(tickLabels)
  ax2.set_xlabel('Wave vector, q (rad / unit distance)',fontsize=16)
  ax2.grid(b=True, which='major', color='k', linestyle='-')

  
plt.tight_layout()
plt.subplots_adjust(hspace=.05)
plt.savefig("../plots/qDependance/V{:d}_T{}_FreePeak_TransDen.svg".format(maxVel,trackLength))
plt.savefig("../plots/qDependance/V{:d}_T{}_FreePeak_TransDen.pdf".format(maxVel,trackLength))
#plt.savefig("../plots/qDependance/V{:d}_d{:.03f}_JamPeak.png".format(maxVel,density))
plt.show()
#raw_input()


