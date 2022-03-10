import numpy as np
import sys
import matplotlib.pyplot as plt

xOffset = float(sys.argv[1])
yOffset = float(sys.argv[2])

thetas = [131.72,132.51,133.30,134.09,134.87,135.66,136.45,137.24,138.03,138.82,139.60,140.39,141.18,141.97,142.76,143.55,144.34,145.12,145.91,146.70,147.49,148.28,149.07,149.85,150.64,151.43,152.22,153.01,153.80,154.59,155.37,156.16,156.95,157.74,158.53,159.32,160.10,160.89,161.68,162.47]
halfthetas = []
for i in range(len(thetas)-1):
  halfthetas.append(thetas[i])
  halfthetas.append((thetas[i]+thetas[i+1])/2)
halfthetas.append(thetas[-1])
offsetRadius = np.sqrt(xOffset**2 + yOffset**2)
print(offsetRadius)

deltaPhi = []
#for theta in thetas:
for theta in halfthetas:
  perp = -3.2*np.tan(theta*np.pi/180)
  minRadius = perp - offsetRadius
  maxRadius = perp + offsetRadius
  if(minRadius >= 1.1 and maxRadius <= 3.5):
    phiRange = [-180.0,180.0]
  elif(minRadius < 1.1):
    phi = -np.pi
    jr = np.sqrt(perp**2 + xOffset**2 + yOffset**2 + 2*perp*xOffset*np.cos(phi) + 2*perp*yOffset*np.sin(phi))
    if(jr < 1.1):
      greater = False
      phiRange = []
    else:
      greater = True
      phiRange = [-180.0]
    for tenphi in range(1,3600):
      phi = (tenphi/10 - 180)*np.pi/180
      jr = np.sqrt(perp**2 + xOffset**2 + yOffset**2 + 2*perp*xOffset*np.cos(phi) + 2*perp*yOffset*np.sin(phi))
      #if(theta > 161.9):
        #print(jr)
      if(jr >= 1.1 and greater == False):
        greater = True
        phiRange.append(phi*180/np.pi)
      elif(jr < 1.1 and greater == True):
        greater = False
        phiRange.append(phi*180/np.pi)
    if phiRange[0] == -180.0:
      phiRange.append(180.0)
  elif(maxRadius > 3.5):
    phi = -np.pi
    jr = np.sqrt(perp**2 + xOffset**2 + yOffset**2 + 2*perp*xOffset*np.cos(phi) + 2*perp*yOffset*np.sin(phi))
    if(jr > 3.5):
      greater = True
      phiRange = []
    else:
      greater = False
      phiRange = [-180.0]
    for tenphi in range(1,3600):
      phi = (tenphi/10 - 180)*np.pi/180
      jr = np.sqrt(perp**2 + xOffset**2 + yOffset**2 + 2*perp*xOffset*np.cos(phi) + 2*perp*yOffset*np.sin(phi))
      if(jr > 3.5 and greater == False):
        greater = True
        phiRange.append(phi*180/np.pi)
      elif(jr <= 3.5 and greater == True):
        greater = False
        phiRange.append(phi*180/np.pi)
    if phiRange[0] == -180.0:
      phiRange.append(180.0)
  if(theta in thetas):
    print(theta)
    print(phiRange)
  if(len(phiRange) == 2):
    deltaPhi.append(phiRange[1] - phiRange[0])
  elif(len(phiRange) == 4):
    deltaPhi.append(phiRange[1] - phiRange[0] + phiRange[3] - phiRange[2])
  else:
    print("No")
print(deltaPhi)
print(len(deltaPhi))
