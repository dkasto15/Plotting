import os;
from numpy import loadtxt
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cmx
import matplotlib.colors as colors
from sys import exit

# # # INSTRUCTIONS: Create a directory in the FocusFinder directory with the name of the investigated sample. Put the TempRamp and PowRamp directories for the sample in the sample directory. Fill in user input below and then run the program. 

# # # User input

path='/home/dkasto/Desktop/FocusFinder' # path to directory with the sample directory and FocusFinder.py
delimit=' ' #delimiter in data files, usually ',', ' ' or '\t'

# # # Find the corresponding sample to a data point
def whichSampleIndex(plPos):
    c,d=0,0
    while plPos > c-1: #because index and not length
        c=c+a[d]
        d=d+1
    return(d-1)
    
# # # Convert the "voltages" to distances from focal point
def distConverter(distArray):
    for i in distArray:
        if i.endswith(".10") or i.endswith(".60"):
            i_corr=-float(i[:3])+0.1
            distFinal.append(i_corr)
        else:
            distFinal.append(float(i))

# # # Fit quality check 
def qual_check(fitCheck):
    errData='********* Failed fits scan *********\n'
    fitCheck=fitCheck[0]
    FWHM=fitCheck[:,3]
    distOrder=fitCheck[:,8]
    distOrder = [str(round(float(i),1))+'0' for i in distOrder]
    for i in distOrder:
        if i.endswith(".10") or i.endswith(".60"):
            i_corr=-float(i[:3])+0.1
            distFinalErr.append(i_corr)
        else:
            distFinalErr.append(float(i))
    samp=samples[0]
    for k in range(len(FWHM)):
        if k==(len(FWHM)-1):
            if FWHM[k]>2*FWHM[k-1]:
                errData=errData + ('Check TempRamp'+str(distOrder[k])+', T='+str(temp)+' K for sample '+str(samp)+'.\n')
    #             # errSum=errSum+1
                errDistNr.append(distOrder[k])
                errEmiss.append(fitCheck[k,1])
            else:
                t=2 #TODO
                # pl_fitCorr_i[k,:]=pl_i[k][:]
        else:
            if FWHM[k]>2*FWHM[k+1]:
                errData=errData + ('Check TempRamp'+str(distOrder[k])+', T='+str(temp)+' K for sample '+str(samp)+'.\n')
    #             # errSum=errSum+1
                errDistNr.append(distOrder[k])
                errEmiss.append(fitCheck[k,1])
            else:
                t=2 #TODO
                # pl_fitCorr_i[k,:]=pl_i[k][:]
    # sum=0
    # pl_fitCorr_i_final=pl_fitCorr_i
    # for k in range(len(FWHM_fit)):
    #     if  pl_fitCorr_i[k,0] == 0:
    #         pl_fitCorr_i_final=np.delete(pl_fitCorr_i_final, (sum), axis=0)
    #     else:
    #         sum=sum+1
    # pl_fitCorr_i_final=pl_fitCorr_i_final.tolist()
    # pl_fitCorr.append(pl_fitCorr_i_final)
    return(errData)
    

# # # Plotting (except for Arrhenius/Tc)

def plotTime(y_col,title,xlabel,ylabel,legendLoc):
    global s
    plt.figure(s)
    s=s+1
    samp=samples[whichSampleIndex(j)]
    title_final=title + 'Sample: ' + str(samp)
    plt.title(title_final)
    plt.plot(distFinalSorted, emissSorted, '.-', label='Sample: '+str(samp))
    plt.plot(errDistNr, errEmiss,'o', markersize=10, fillstyle='none', markeredgewidth=2, color='red')
    plt.legend(loc=legendLoc)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)


# # # Loading pl files
os.chdir(path) 

theFile='./PLmaxsig.txt'
theFile3='./3dData.txt'

dist,pl,samples,a,diffDist,fitCheck=[],[],[],[],[],[]
for file in os.listdir(path):
    if not file.endswith(".py"):
        samples.append(file)
        sum=0
        os.chdir(path + '/' + file)
        for file2 in os.listdir(os.getcwd()):
            if file2.startswith("Temp"):
                sum=sum+1
                dist.append(file2[-5:-1])
                os.chdir(os.getcwd()+'/TempRamp'+file2[-5:-1]+'V')
                pl.append(loadtxt(theFile, comments="#", delimiter=delimit, unpack=False))
                CO2_check.append(loadtxt(theFile3, comments="#", delimiter=delimit, unpack=False))
                os.chdir(path + '/' + file)
            if file2.startswith("PowRamp"):
                #fitFail.append(file2[-5:-1])
                os.chdir(os.getcwd()+'/PowRamp'+file2[-4:-1]+'K')
                fitCheck.append(loadtxt(theFile, comments="#", delimiter=delimit, unpack=False))
                os.chdir(path + '/' + file)
        a.append(sum) #num of voltages for each sample

# # Sorting etc.

distFinal=[]
distConverter(dist)
temp=int(pl[0][0])

emiss=[]
for j in range(len(pl)):
    emiss.append(pl[j][1])

distFinalSorted=sorted(distFinal)
emissSorted=np.ones(len(distFinalSorted))
for k in distFinalSorted:
    emissSorted[distFinalSorted.index(k)]=emiss[distFinal.index(k)]

# # Fit check
distFinalErr=[]
errDistNr=[]
errEmiss=[]
print(qual_check(fitCheck))
# fitCheck=fitCheck[0]
# distOrder=fitCheck[:,8]
# distOrder = [str(round(float(i),1))+'0' for i in distOrder]
# for i in distOrder:
#     if i.endswith(".10") or i.endswith(".60"):
#         i_corr=-float(i[:3])+0.1
#         distFinalErr.append(i_corr)
#     else:
#         distFinalErr.append(float(i))


# errSum=0
samp=samples[0]

# # Plotting

s=1 #figure number
plotTime(1,'Focal point scan - ','Position elec. stage [mm]','PL emission \n Wavenumber [cm$^{-1}$]',4)
    
plt.show()