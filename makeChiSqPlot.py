"""This program reads the .chisq output of the chisqsurface program and makes plots of the results. It will generate the following four plots in the same directory as the .chisq file: target chisq surface, projectile chisq surface, total chisq surface, and total chisq surface with a one-sigma cut. This program takes two inputs, the .chisq file (required) and the energy levels for each matrix element (optional). To use the optional input, pass the --levels flag followed by the four levels in the format xpy, where x is the spin, p is the parity (+ or -) and y is the number of the state (i.e. first 2+ would be 2+1, second 3- would be 3-2, etc.). Each level should be seperated by a space. The first two levels will be used to make the axis label for the x axis, and the last two will be used for the y axis. If you do not supply the levels, the default axis labels will be <2+1|E2|0+1> for the x axis and <2+1|E2|2+1> for the y axis."""
import numpy as np
import argparse
import matplotlib.pyplot as plt

#define the command line arguments for this program.
parser = argparse.ArgumentParser(description="Plot chisq surfaces")
parser.add_argument("chisqFile",type=str,nargs='?',help='The input file with the chisq values')
parser.add_argument("--levels",type=str,nargs=4,default=[2+1,0+1,2+1,2+1],help='The transition states in question, fomatted as spin parity number, so 2+1 is the first 2+ state')

#read the arguments and open the chisq file
args = parser.parse_args()
chisqFile = open(args.chisqFile,'r')

#Each line in the chisq file has five tab-separated entries: x matrix element, y matrix element, projectile chisq, target chisq, and total chisq. We will read them into lists here.
me1 = []
me2 = []
projChi = []
targChi = []
totChi = []

line = chisqFile.readline()
while line:
  splitLine = line.split("\t")
  me1.append(float(splitLine[0]))
  me2.append(float(splitLine[1]))
  projChi.append(float(splitLine[2]))
  targChi.append(float(splitLine[3]))
  totChi.append(float(splitLine[4][:-1]))
  line = chisqFile.readline()

#In order to make the 2D plots, we need to find the number of "unique" entries in each list of matrix elements
unique1 = []
unique2 = []
for i in range(len(me1)):
  if me1[i] not in unique1:
    unique1.append(me1[i])
  if me2[i] not in unique2:
    unique2.append(me2[i])

#Sort the matrix elements. The y axis is sorted in reverse order for plotting purposes.
unique1.sort()
unique2.sort(reverse=True)

#Now we make four matrices to store the chisq values. These matrices will be the numerical representation of our 2D plots. We also want to track the minimum chisq value for the one sigma surface, and we also want to track the bounds on the matrix elements in that one sigma surface so we can extract the uncertainties.
projChiMatrix = np.zeros((len(unique2),len(unique1)))
targChiMatrix = np.zeros((len(unique2),len(unique1)))
totChiMatrix = np.zeros((len(unique2),len(unique1)))
chiPlusOneMatrix = np.zeros((len(unique2),len(unique1)))

chiMin = min(totChi)
me1min = 999
me1max = -999
me2min = 999
me2max = -999

for i in range(len(me1)):
  xIndex = unique1.index(me1[i])
  yIndex = unique2.index(me2[i])
  projChiMatrix[yIndex,xIndex] = projChi[i]
  targChiMatrix[yIndex,xIndex] = targChi[i]
  totChiMatrix[yIndex,xIndex] = totChi[i]
#If the totcal chisq is within the one sigma surface, we check and update the error bounds
  if(totChi[i] <= chiMin + 1):
    chiPlusOneMatrix[yIndex,xIndex] = totChi[i]
    if(me1[i] < me1min):
      me1min = me1[i]
    elif(me1[i] > me1max):
      me1max = me1[i]
    if(me2[i] < me2min):
      me2min = me2[i]
    elif(me2[i] > me2max):
      me2max = me2[i]
    if(totChi[i] == chiMin):
      me1best = me1[i]
      me2best = me2[i]

#Plots
extent = np.min(unique1), np.max(unique1), np.min(unique2), np.max(unique2)

plt.imshow(projChiMatrix,cmap='hot',interpolation='none',extent=extent,aspect='auto')
#Use the input arguments for the levels to generate the axis labels.
plt.xlabel("<$%s^%s_%s$|E2|$%s^%s_%s$> eb" % (args.levels[0][0], args.levels[0][1], args.levels[0][2], args.levels[1][0], args.levels[1][1], args.levels[1][2]))
plt.ylabel("<$%s^%s_%s$|E2|$%s^%s_%s$> eb" % (args.levels[2][0], args.levels[2][1], args.levels[2][2], args.levels[3][0], args.levels[3][1], args.levels[3][2]))
plt.colorbar()
plt.title("Projectile $\u03A7^2$")
plt.savefig(args.chisqFile.split(".")[0] + str("_proj.png"))
plt.close()

plt.imshow(targChiMatrix,cmap='hot',interpolation='none',extent=extent,aspect='auto')
plt.xlabel("<$%s^%s_%s$|E2|$%s^%s_%s$> eb" % (args.levels[0][0], args.levels[0][1], args.levels[0][2], args.levels[1][0], args.levels[1][1], args.levels[1][2]))
plt.ylabel("<$%s^%s_%s$|E2|$%s^%s_%s$> eb" % (args.levels[2][0], args.levels[2][1], args.levels[2][2], args.levels[3][0], args.levels[3][1], args.levels[3][2]))
plt.colorbar()
plt.title("Target $\u03A7^2$")
plt.savefig(args.chisqFile.split(".")[0] + str("_targ.png"))
plt.close()

plt.imshow(totChiMatrix,cmap='hot',interpolation='none',extent=extent,aspect='auto')
plt.xlabel("<$%s^%s_%s$|E2|$%s^%s_%s$> eb" % (args.levels[0][0], args.levels[0][1], args.levels[0][2], args.levels[1][0], args.levels[1][1], args.levels[1][2]))
plt.ylabel("<$%s^%s_%s$|E2|$%s^%s_%s$> eb" % (args.levels[2][0], args.levels[2][1], args.levels[2][2], args.levels[3][0], args.levels[3][1], args.levels[3][2]))
plt.colorbar()
plt.title("Total $\u03A7^2$")
plt.savefig(args.chisqFile.split(".")[0] + str("_tot.png"))
plt.close()

data_masked = np.ma.masked_where(chiPlusOneMatrix==0,chiPlusOneMatrix)

plt.imshow(data_masked,cmap='hot',interpolation='none',extent=extent,aspect='auto')
plt.xlabel("<$%s^%s_%s$|E2|$%s^%s_%s$> eb" % (args.levels[0][0], args.levels[0][1], args.levels[0][2], args.levels[1][0], args.levels[1][1], args.levels[1][2]))
plt.ylabel("<$%s^%s_%s$|E2|$%s^%s_%s$> eb" % (args.levels[2][0], args.levels[2][1], args.levels[2][2], args.levels[3][0], args.levels[3][1], args.levels[3][2]))
plt.colorbar()
plt.title("$\u03A7^2$+1 Cut")
plt.savefig(args.chisqFile.split(".")[0] + str("_tot_onesigma.png"))
plt.close()

#print the bounds of the chisq+1 surface
print("Matrix Element One Best Fit: (" + str(me1best) + ")")
print("Matrix Element One Bounds: (" + str(me1min) + "," + str(me1max) + ")")
print("Matrix Element Two Best Fit: (" + str(me2best) + ")")
print("Matrix Element Two Bounds: (" + str(me2min) + "," + str(me2max) + ")")
