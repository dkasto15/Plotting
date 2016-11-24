import os;
from numpy import loadtxt
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cmx
import matplotlib.colors as colors
from sys import exit

# # # INSTRUCTIONS: Create directories in the SampleComparisons directory with the names of the samples. Put the TempRamp directories for each sample in the corresponding sample directory. Fill in user input below and then run the program. 

# # # User input

path='/home/dkasto/Desktop/SampleComparisons_v2' # path to directory with sample directories and SampleComparisons.py
distFromFoc=3 #Distance from the focal point for the pump laser (in mm)

integPL='y' #Write 'y' if you want integrated PL graphs
emissPL='y' #Write 'y' if you want PL emission graphs
FWHM_PL='y' #Write 'y' if you want FWHM PL emission graphs
Tc='' #Write 'y' if you want Tc graphs
errFitCorr='y' #Write 'y' if you want the graphs to mark data points with erroneous Voigt fits
CO2Alert='y' #Write 'y' if you want the graphs to mark data points with Voigt fits close to the CO2 absorption region
renameSamples='' # Present samples with e.g. a specific property instead of their sample names in the legend of the plots. Example: "5 QWs" instead of "Sample 1718"


# # # Different colors for the different curves

def get_cmap(N):
    '''Returns a function that maps each index in 0, 1, ... N-1 to a distinct 
    RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv') 
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color

# # # Find the corresponding sample to a data point
def whichSampleIndex(plPos):
    c,d=0,0
    while plPos > c-1: #because index and not length
        c=c+a[d]
        d=d+1
    return(d-1)
    
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
    
# # # Rename a sample in the legend of the plots
def renameASample(sample):
    b=samples.index(sample)
    return(renamedSamp[b])

# # # Fit quality check 

def qual_check(pl,a):
    errData='********* Failed fits scan *********\n'
    for i in range(len(pl)):
        pl_i=np.asarray(pl[i])
        pl_fitCorr_i=np.zeros((len(pl_i[:,0]),len(pl_i[0,:])))
        FWHM_fit=[]
        gamma_fit=[]
        # errSum=0
        samp=samples[whichSampleIndex(i)]
        for j in range(len(pl_i[:,0])):
            FWHM_fit.append(pl_i[j,3])
            gamma_fit.append(pl_i[j,5])
        for k in range(len(FWHM_fit)):
            if k==(len(FWHM_fit)-1):
                if FWHM_fit[k]>2*FWHM_fit[k-1] or abs(gamma_fit[k])>abs(3*gamma_fit[k-1]):
                    errData=errData + ('Check TempRamp'+str(volt[i])+', T='+str(pl[i][k,0])+' K for sample '+str(samp)+'.\n')
                    # errSum=errSum+1
                    errTempNr.append(k)
                    errPlNr.append(i)
                else:
                    pl_fitCorr_i[k,:]=pl_i[k][:]
            else:
                if FWHM_fit[k]>2*FWHM_fit[k+1] or abs(gamma_fit[k])>abs(3*gamma_fit[k+1]):
                    errData=errData+('Check TempRamp'+str(volt[i])+', T='+str(pl[i][k,0])+' K for sample '+str(samp)+'.\n')
                    # errSum=errSum+1
                    errTempNr.append(k)
                    errPlNr.append(i)
                else:
                    pl_fitCorr_i[k,:]=pl_i[k][:]
        sum=0
        pl_fitCorr_i_final=pl_fitCorr_i
        for k in range(len(FWHM_fit)):
            if  pl_fitCorr_i[k,0] == 0:
                pl_fitCorr_i_final=np.delete(pl_fitCorr_i_final, (sum), axis=0)
            else:
                sum=sum+1
        pl_fitCorr_i_final=pl_fitCorr_i_final.tolist()
        pl_fitCorr.append(pl_fitCorr_i_final)
    return(errData)
    
# # # CO2 area alert

def CO2_alert():
    alertMessage='********* CO2 area alert *********\n'
    for i in range(len(CO2_check)):
        for n in range(len(pl[0][:,0])):
            j=pl[0][n,0] #temp instead of index
            f = sorted(np.where(CO2_check[i][:,0] == j))
            g=CO2_check[i][f[0],1:3] 
            h=np.argpartition(g[:,1], -5)[-5:] #yields the indices to the 5 biggest values (intensities) in g
            h=g[h,0] #corresponding wavenumbers
            inRange=np.asarray([int(l) for l in [k > 2250 and k < 2450  for k in h]])
            elemInRange=0
            for m in inRange:
                elemInRange=elemInRange+m
            if elemInRange > 3:
                samp=samples[whichSampleIndex(i)]
                alertMessage=alertMessage+'Sample '+str(samp)+': Temp='+str(j)+', Volt='+str(volt[i])+' might be a faulty bastard.\n'
                badFit=0
                if len(errTempNr)>0:
                    for b in range(len(errTempNr)): #Append values if the data point isn't from a failed fit
                        if errTempNr[b]==n and errPlNr[b]==i:
                            badFit=badFit+1
                    if badFit == 0:
                        errTempNrCO2.append(n)
                        errPlNrCO2.append(i)                            
                else:
                    errTempNrCO2.append(n)
                    errPlNrCO2.append(i)
    return(alertMessage)

# # # Plotting (except for Arrhenius/Tc)

def plotTime(y_col,title,xlabel,ylabel,legendLoc):
    global s
    for i in diffPowers:
        plt.figure(s)
        s=s+1
        title_final=title + str(i) +' V\nDistance from focal point='+str(distFromFoc)+' mm'
        plt.title(title_final)
        b = np.where(volt == i)[0]
        colSum=0
        for j in b:
            samp=samples[whichSampleIndex(j)]
            if renameSamples=='y':
                sampLabel=renameASample(samp)
            else:
                sampLabel='Sample: '+str(samp)
            col = cmap(colSum+1)
            colSum=colSum+1
            plt.plot(pl[j][:,0], pl[j][:,y_col], '.-', label=sampLabel, color=col)
            plt.legend(loc=legendLoc)
            plt.ylabel(ylabel)
            plt.xlabel(xlabel)

def plotTimeErrCorr(y_col,title,xlabel,ylabel,legendLoc): #Also includes CO2
    global s
    global CO2Alert
    for i in diffPowers:
        plt.figure(s)
        s=s+1
        title_final=title + str(i) +' V\nDistance from focal point='+str(distFromFoc)+' mm'
        plt.title(title_final)
        b = np.where(volt == i)[0]
        colSum=0
        for j in b:
            samp=samples[whichSampleIndex(j)]
            if renameSamples=='y':
                sampLabel=renameASample(samp)
            else:
                sampLabel='Sample: '+str(samp)
            col = cmap(colSum+1)
            colSum=colSum+1
            pl_fitCorr_j=np.array(pl_fitCorr[j])
            plt.plot(pl_fitCorr_j[:,0], pl_fitCorr_j[:,y_col], '.-', label=sampLabel, color=col)
            e=np.where(errPlNr == j)[0]
            if np.size(e) > 0:  #Does this sample+voltage contain at least one critical data point
                for k in e:
                    if errTempNr[k] > 0:
                        plt.plot(pl[j][errTempNr[k],0], (pl[j][errTempNr[k]-1,y_col]+pl[j][errTempNr[k]+1,y_col])/2, 'o', markersize=10, label='Deleted data point (failed fit)', color=col)
                    else:
                        plt.plot(pl[j][errTempNr[k],0], pl[j][errTempNr[k]+1,y_col], 'o', markersize=10, label='Deleted data point (failed fit)', color=col)
            if CO2Alert=='y':
                e=np.where(errPlNrCO2 == j)[0]
                if np.size(e) > 0: #Does this sample+voltage contain a critical data point
                    labelOnlyOnce=0
                    for k in e:
                        if labelOnlyOnce==0: #Only one label per colour (sample)
                            labelOnlyOnce=labelOnlyOnce+1
                            if errTempNrCO2[k] > 0:
                                plt.plot(pl[j][errTempNrCO2[k],0], pl[j][errTempNrCO2[k],y_col], 'o', markersize=10, fillstyle='none', markeredgewidth=2, label='Critical data point (CO2)', color=col)
                            else:
                                plt.plot(pl[j][errTempNrCO2[k],0], pl[j][errTempNrCO2[k],y_col], 'o', markersize=10, fillstyle='none', markeredgewidth=2, label='Critical data point (CO2)', color=col)
                        else:
                            if errTempNrCO2[k] > 0:
                                plt.plot(pl[j][errTempNrCO2[k],0], pl[j][errTempNrCO2[k],y_col],'o', markersize=10, fillstyle='none', markeredgewidth=2, color=col)
                            else:
                                plt.plot(pl[j][errTempNrCO2[k],0], pl[j][errTempNrCO2[k],y_col],'o', markersize=10, fillstyle='none', markeredgewidth=2, color=col)
            plt.legend(loc=legendLoc)
            plt.ylabel(ylabel)
            plt.xlabel(xlabel)

# # # Loading pl files
os.chdir(path) 

theFile='./PLmaxsig.txt'
theFile2='./ArrValues.txt'
theFile3='./3dData.txt'

numOfSamp=0
volt,pl,samples,a,diffPowers,arrh,samplesArrh,CO2_check=[],[],[],[],[],[],[],[]
for file in os.listdir(path):
    if not file.endswith(".py"):
        numOfSamp=numOfSamp+1
        sum=0
        samples.append(file)
        os.chdir(path + '/' + file)
        for file2 in os.listdir(os.getcwd()):
            if file2.startswith("Temp"):
                sum=sum+1
                volt.append(file2[-5:-1])
                os.chdir(os.getcwd()+'/TempRamp'+file2[-5:-1]+'V')
                delimit=whichDelimiter(theFile)
                pl.append(loadtxt(theFile, comments="#", delimiter=delimit, unpack=False))
                delimit=whichDelimiter(theFile3)
                CO2_check.append(loadtxt(theFile3, comments="#", delimiter=delimit, unpack=False))
                if file2[-5:-1] not in diffPowers:
                    diffPowers.append(file2[-5:-1])
                os.chdir(path + '/' + file)
            else:
                if file2.startswith("Arr"):
                    os.chdir(os.getcwd()+'/ArrPlot')
                    samplesArrh.append(file)
                    delimit=whichDelimiter(theFile2)
                    arrh.append(loadtxt(theFile2, comments="#", delimiter=delimit, unpack=False))
                os.chdir(path + '/' + file)
        a.append(sum) #num of voltages for each sample

# # #

cmap = get_cmap(len(pl)+1)

diffPowers=sorted(diffPowers)
diffPowers = [float(i) for i in diffPowers]
volt = [float(i) for i in volt]
volt=np.asarray(volt)
s=1 #figure number

# # # Fit quality check and CO2 alert
pl_fitCorr=[] #pl without data points from failed fits
#errVolt=[]
errTempNr=[]
errPlNr=[]
print(qual_check(pl,a))
errTempNr=np.asarray(errTempNr)
errPlNr=np.asarray(errPlNr)

#errVoltCO2=[]
errTempNrCO2=[]
errPlNrCO2=[]
print(CO2_alert())
errTempNrCO2=np.asarray(errTempNrCO2)
errPlNrCO2=np.asarray(errPlNrCO2)

# # # Rename samples in legend of plots
if renameSamples=='y':
    renamedSamp=[]
    for i in samples:
        renamedSamp.append(input('What would you like to call Sample '+str(i)+' in the legend of the plots?\n'+'--> '))
        


# # #

if errFitCorr=='y':
    
    # # integPL errFitCorr
    if integPL=='y':
        plotTimeErrCorr(2,'Corrected Integrated PL - ','Temperature [K]','Integrated PL intensity [a.u]',1)
    
    # # emissPL errFitCorr
    if emissPL=='y':
        plotTimeErrCorr(1,'Corrected PL emission - ','Temperature [K]','PL emission \n Wavenumber [cm$^{-1}$]',4)
                
    # # FWHM_PL errFitCorr
    if FWHM_PL=='y':
        plotTimeErrCorr(3,'Corrected PL emission FWHM - ','Temperature [K]','PL emission FWHM \n Wavenumber [cm$^{-1}$]',4)
    
    # # Tc
    if Tc=='y':
        print('***\n N.B. To get an adjusted Tc graph, remove the erroneous measurements pointed out above from the PLmaxsig files before you run the ArrPlot program, then run this program without errFitCorr.\n***\n')

else:
    # # integPL
    if integPL=='y':
        plotTime(2,'Integrated PL - ','Temperature [K]','Integrated PL intensity [a.u]',1)
    
    # # emissPL
    if emissPL=='y':
        plotTime(1,'PL emission - ','Temperature [K]','PL emission \n Wavenumber [cm$^{-1}$]',4)
    
    # # FWHM_PL
    if FWHM_PL=='y':
        plotTime(3,'PL emission FWHM - ','Temperature [K]','PL emission FWHM \n Wavenumber [cm$^{-1}$]',4)
        
    # # Tc
    if Tc=='y':
        plt.figure(s)
        plt.title('Tc for different powers and samples')
        colSum=0
        incr=0
        for k in samplesArrh:
            col = cmap(colSum+1)
            colSum=colSum+1
            samp=k
            if renameSamples=='y':
                sampLabel=renameASample(samp)
            else:
                sampLabel='Sample: '+str(samp)
            plt.plot(arrh[incr][:,0], arrh[incr][:,5], '.-', label=sampLabel, color=col)
            incr=incr+1
            plt.legend()
            plt.ylabel('Critical temperature Tc [K]')
            plt.xlabel('Excitation voltage [V]')


plt.show()