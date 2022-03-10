"""The purpose of this code is to generate a database of the expected sector distributions for various beam spots. The distributions can be compared to the distribution measured during an experiment to extract the beam spot. The distributions are generated via Monte Carlo simulation."""
import numpy as np
import bisect
import matplotlib.pyplot as plt

#The phiBounds array stores the angular coverage of each sector, starting from 0 degrees. the sectorMap array maps these angle ranges to the actual sector numbers. thetaBounds aand ringMap do the same thing for the rings.
phiInc = np.pi/16
phiBounds = [x*phiInc for x in range(33)]
sectorMap = [9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,1,2,3,4,5,6,7,8]
thetaBounds = [np.arctan((1.1+0.1*x)/2.8) for x in range(25)]
ringMap = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]

#Open the database file that we will write to
database = open("beamspotDatabase.txt",'w',buffering=1)

#Alpha and beta are coefficients used for integrating the Rutherford cross section for simulating thetas. alpha and beta are equal to the csc^2(theta/2) for the minimum and maximum theta angles we want to simulate, respectively. Since we are looking at the downstream detector, we can set the maximum angle to 90 degrees. However, because the cross section diverges at zero degrees, the minimum integration angle is set to exactly the smallest angle we could ever measure when accounting for the array geometry and the largest beamspot we are simulating. A smaller angle could also be chosen, but due to the diverging nature of the cross section, this will result in a large number of MC points that don't hit the detector.
alpha = 66.9848
beta = 5.62012

#The offsets we look at here cover the range from -2mm to 2mm, in increments of 0.1 mm
xOffsets = [x/100. - 0.2 for x in range(41)]
yOffsets = [x/100. - 0.2 for x in range(41)]
residuals = np.zeros((len(xOffsets),len(yOffsets)))
iters = 200000

#Loop over all combinations of x and y
for x in xOffsets:
  for y in yOffsets:
#store the sector and ring distributions
    mcSectorDist = [0]*32
    mcRingDist = [0]*24
    #print(x,y)
    for j in range(iters):
#Randomly generate a theta by integrating the Rutherford differential cross section
      rand_theta = np.random.random()
      mc_theta = 2*np.arcsin((alpha + (beta-alpha)*rand_theta)**(-0.5))
#Randomly generate a phi from a uniform distribution
      mc_phi = (np.random.random()*2-1)*np.pi
#Rho is the distance from the beam axis to the point of interest. This is determined only by theta. The 2.8 used here is the 2.8 cm from the target to the downstream detector.
      mc_rho = 2.8*np.tan(mc_theta)
#Using rho and phi, we can compute x and y, and then apply the offsets to get the x,y coordinates in the JANUS frame
      mc_x = mc_rho*np.cos(mc_phi)
      mc_y = mc_rho*np.sin(mc_phi)
      off_x = mc_x - x
      off_y = mc_y - y
#Now we convert the JANUS coordinates from Cartesian to rho, theta, and phi coordinates in the plane of the detector
      off_rho = np.sqrt(off_x**2 + off_y**2)
      off_theta = np.arctan2(off_rho,2.8)
      off_phi = np.arctan2(off_y,off_x)
      if off_phi < 0:
        off_phi += 2*np.pi
#We use the bisect function and our sector and ring maps to identify the sector and ring where the particle would hit
      index = bisect.bisect_left(thetaBounds,off_theta)-1
      if(index >= 0 and index < 24):
        ring = ringMap[index]
        mcRingDist[ring-1] += 1
        index = bisect.bisect_left(phiBounds,off_phi)-1
        sector = sectorMap[index]
        mcSectorDist[sector-1] += 1
#Normalize the ring and sector distributions and then right them to the database
    validCounts = sum(mcSectorDist)
    mcSectorDist = [x/validCounts for x in mcSectorDist]
    mcRingDist = [x/validCounts for x in mcRingDist]
    database.write(str(x)+","+str(y)+"\n")
    sDistWrite = ""
    for p in mcSectorDist:
      sDistWrite += str(p) + ","
    sDistWrite = sDistWrite[:-1] + "\n"
    database.write(sDistWrite)
    rDistWrite = ""
    for p in mcRingDist:
      rDistWrite += str(p) + ","
    rDistWrite = rDistWrite[:-1] + "\n"
    database.write(rDistWrite)

database.close()
