"""This function generates the phi patter needed to apply the beam spot corrections to GOSIA2. It takes as input (in order): the x offset (mm), the y offset (mm), the detector number (1 for downstream, 0 for upstream), and the number of theta increments that will be used in GOSIA2's OP,INTI. The function will return the theta increments, the list of phi ranges for each theta, and the total phi coverage at each value of theta and the intermediate value between each pair of thetas."""
import numpy as np
import sys
import matplotlib.pyplot as plt
import argparse

#Take in all of the arguments
parser = argparse.ArgumentParser(description="Generate Inputs for GOSIA2 beam spot corrections")
parser.add_argument("xOffset",type=float,nargs='?',help='X offset for the beamspot')
parser.add_argument("yOffset",type=float,nargs='?',help='Y offset for the beamspot')
parser.add_argument("detnum",type=int,nargs='?',help='JANUS detector number')
parser.add_argument("nInc",type=int,nargs='?',help='Number of theta increments desired')

args = parser.parse_args()
xOffset = args.xOffset
yOffset = args.yOffset
detnum = args.detnum
nInc = args.nInc

#For the downstream detector, the z offset is 2.8 mm. For the upstream detector, it's 3.2 mm. 
if detnum == 1:
  zAxis = 2.8
elif detnum == 0:
  zAxis = -3.2
else:
  print("ERROR: Invalid detector number. Detector number must be 1 for downstream detector and 0 for upstream detector.")

#Compute the theta values at which we will find the phi cooverage.
#Rho is the total distance in the x-y plane. rhoMin and rhoMax are the smallest and largest possible values or rho that will hit the detector when accounting for the offsets.
rhoOffset = np.sqrt(xOffset**2 + yOffset**2)
rhoMin = 1.1 - rhoOffset
rhoMax = 3.5 + rhoOffset

#Using rhoMin and rhoMax, we ccmpute the minimum and maximum theta angles covered by the detector. We then use this plus the number of increments specified by the user to generate the list of theta values. We apply a small offset to thetaMin and thetaMax to prevent precision errors.
if detnum == 1:
  thetaMin = np.arctan2(rhoMin,zAxis)*180.0/np.pi + 0.05
  thetaMax = np.arctan2(rhoMax,zAxis)*180.0/np.pi - 0.05
elif detnum == 0:
  thetaMin = np.arctan2(rhoMax,zAxis)*180.0/np.pi + 0.05
  thetaMax = np.arctan2(rhoMin,zAxis)*180.0/np.pi - 0.05

thetaInc = (thetaMax - thetaMin)/(nInc - 1)
thetas = [thetaMin + thetaInc*x for x in range(nInc)]
print(*["%.2f" % x for x in thetas], sep=',')

#The total phi coverage needs to be computed at not only the values in thetas, but also each intermediate value between theta increments. We use the halfthetas array for this.
halfthetas = []
for i in range(len(thetas)-1):
  halfthetas.append(thetas[i])
  halfthetas.append((thetas[i]+thetas[i+1])/2)
halfthetas.append(thetas[-1])

#deltaPhi stores the total angular coverage for each value in halfthetas
deltaPhi = []
for theta in halfthetas:
#rho is the total distance in the x-y plane with respect to the beam axis. We then apply the offsets to find the minimum and maximum rho in the detector frame.
  rho = zAxis*np.tan(theta*np.pi/180)
  minRadius = rho - rhoOffset
  maxRadius = rho + rhoOffset
#If the min and max rho are always within the coverage of the detecotr, the phi range is all phi and we are done.
  if(minRadius >= 1.1 and maxRadius <= 3.5):
    phiRange = [-180.0,180.0]
#If the minimum rho is less than the detector minimum, we need to compute the phi coverage. We will start at -180 degrees and iterate over phi in increments of 0.1 degrees to find all of the points where rho crosses the edge of the detector. This will give us the total phi coverage. If -180 degrees is covered by the detector, we must count this as a crossing point and then also count +180 degrees as a crossing point, because that's how it has to be input into GOSIA2. 
  elif(minRadius < 1.1):
    phi = -np.pi
    jr = np.sqrt(rho**2 + xOffset**2 + yOffset**2 + 2*rho*xOffset*np.cos(phi) + 2*rho*yOffset*np.sin(phi))
    if(jr < 1.1):
      greater = False
      phiRange = []
    else:
      greater = True
      phiRange = [-180.0]
    for tenphi in range(1,3600):
      phi = (tenphi/10 - 180)*np.pi/180
      jr = np.sqrt(rho**2 + xOffset**2 + yOffset**2 + 2*rho*xOffset*np.cos(phi) + 2*rho*yOffset*np.sin(phi))
      if(jr >= 1.1 and greater == False):
        greater = True
        phiRange.append(phi*180/np.pi)
      elif(jr < 1.1 and greater == True):
        greater = False
        phiRange.append(phi*180/np.pi)
#If -180 is covered, +180 is also covered
    if phiRange[0] == -180.0:
      phiRange.append(180.0)
#This is the same as the last block of code, but for the outer edge of the detector instead of the inner edge.
  elif(maxRadius > 3.5):
    phi = -np.pi
    jr = np.sqrt(rho**2 + xOffset**2 + yOffset**2 + 2*rho*xOffset*np.cos(phi) + 2*rho*yOffset*np.sin(phi))
    if(jr > 3.5):
      greater = True
      phiRange = []
    else:
      greater = False
      phiRange = [-180.0]
    for tenphi in range(1,3600):
      phi = (tenphi/10 - 180)*np.pi/180
      jr = np.sqrt(rho**2 + xOffset**2 + yOffset**2 + 2*rho*xOffset*np.cos(phi) + 2*rho*yOffset*np.sin(phi))
      if(jr > 3.5 and greater == False):
        greater = True
        phiRange.append(phi*180/np.pi)
      elif(jr <= 3.5 and greater == True):
        greater = False
        phiRange.append(phi*180/np.pi)
    if phiRange[0] == -180.0:
      phiRange.append(180.0)
#If theta is in thetas, we print the phi coverage of the detector in the format that it has to be input into GOSIA2
  if(theta in thetas):
    print(int(len(phiRange)/2))
    print(*["%.1f" % x for x in phiRange], sep=',')
#We also need to track the total coverage for each iteration, so we do that here. 
  if(len(phiRange) == 2):
    deltaPhi.append(phiRange[1] - phiRange[0])
  elif(len(phiRange) == 4):
    deltaPhi.append(phiRange[1] - phiRange[0] + phiRange[3] - phiRange[2])
  else:
    print("ERROR: phi range does not have the correct number of entries for theta = " + str(theta))
#outDeltaPhi = ["%f.1" % x for x in deltaPhi]
print(*["%.1f" % x for x in deltaPhi], sep=',')
