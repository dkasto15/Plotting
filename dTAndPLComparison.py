import os;
from numpy import loadtxt
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cmx
import matplotlib.colors as colors
from sys import exit

# # # INSTRUCTIONS: blablalba 

# # # User input

path='/home/dkasto/Desktop/dTAndPLComparison' # path to directory with sample directories and SampleComparisons.py
upNoiseLim=3000 #set region you want to plot [wavenumber]
lowNoiseLim=2385 #set region you want to plot [wavenumber]

# # # Plotting
def plotTime(title,xlabel,ylabel,legendLoc):
    global s
    global upNoiseLim
    global lowNoiseLim
    colSum=0
    for k in range(int(len(dT)/2)):
        plt.figure(s)
        title_final=title# + str(i) +' V'
        plt.title(title_final)
        samp=samples[k]
        sampLabel='Sample: '+str(samp)
        col = cmap(colSum+1)
        colSum=colSum+1
        x=dT[k][:,0]
        if x[1]>x[0]: # ascending or descending order
            i=0
            while x[i]<lowNoiseLim:
                i=i+1
            j=i
            while x[j]<upNoiseLim:
                j=j+1
            x=x[i:j]
            y=dT[2*k+1][i:j,1]-dT[2*k][i:j,1]
        if x[1]<x[0]: # ascending or descending order
            i=0
            while x[i]>upNoiseLim:
                i=i+1
            j=i
            while x[j]>lowNoiseLim:
                j=j+1
            x=x[i:j]
            y=dT[2*k+1][i:j,1]-dT[2*k][i:j,1]
        print('ymax at x='+str(x[np.where(y==max(y))[0][0]]))
        # for i in range(2,len(y)-3): # averaging
        #     y[i]=np.mean(y[i-2:i+3]) # averaging
        plt.plot(x, y, '.-', label=sampLabel, color=col)
        plt.legend(loc=legendLoc)
        plt.ylabel(ylabel+temp[2*k+1]+'-T='+temp[2*k]+')')
        plt.xlabel(xlabel)
        plt.xlim(lowNoiseLim,upNoiseLim)

# # # Upper temperature first for each sample
def upTempFirst():
    for k in range(int(len(dT)/2)):
        a=temp[2*k+1]
        b=temp[2*k]
        c=dT[2*k+1]
        d=dT[2*k]
        if a<b:
            temp[2*k+1]=b
            temp[2*k]=a
            dT[2*k+1]=d
            dT[2*k]=c  

# # # Find the delimiter in the data files
def whichDelimiter(dataFile):
    delim=' '
    delimAlt=[',','\t']
    sum=0
    while True:
        try:
            loadtxt(dataFile, comments="#", delimiter=delim, unpack=False)
            break
        except(ValueError):
            delim=delimAlt[sum]
            sum=sum+1
    return(delim)

# # # Different colors for the different curves
def get_cmap(N):
    '''Returns a function that maps each index in 0, 1, ... N-1 to a distinct 
    RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv') 
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color

# # # Loading dT files

dT=[]
samples=[]
temp=[]
numOfSamp=0
for file in os.listdir(path):
    if not file.endswith(".py"):
        numOfSamp=numOfSamp+1
        sum=0
        samples.append(file)
        os.chdir(path + '/' + file)
        for fileDPT in os.listdir(os.getcwd()):
            if fileDPT.endswith(".DPT"):
                pre, ext = os.path.splitext(fileDPT)
                os.rename(fileDPT, pre + '.txt')
        for file2 in os.listdir(os.getcwd()):
            sum=sum+1
            temp.append(file2[8:file2.index('C')])
            # os.chdir(os.getcwd()+'/TempRamp'+file2[-5:-1]+'V')
            delimit=whichDelimiter(file2)
            dT.append(loadtxt(file2, comments="#", delimiter=delimit, unpack=False))
        os.chdir(path)        

# # # 
upTempFirst()
s=1
cmap = get_cmap(len(dT)+1)
plotTime('Differential transmission','Wavenumber [cm$^{-1}$]','Differential transmission (T=',4)
plt.show()