from numpy import array,std,var,average,linspace,sin,cos,pi,sqrt
import numpy as np
import csv
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
from scipy import optimize
from scipy.special import erf
from matplotlib.ticker import NullFormatter
from math import factorial
import matplotlib.mlab as mlab
from scipy.misc import factorial
from scipy.optimize import curve_fit

print("working")

# Parameters
CallAll = False # prints and plots all lab results
Test = False

plotA = True
printA = True
saveA = False

figS = (30,15) # figure size for plots
fontS = 20
LW = 5
labS = 20 # plot axis label size
titS = 24 # plot title size
binN = 50
N_plot = 1001

if CallAll == True:
    plotA = True
    printA = True
    saveA = True
if Test == True:
    plotA = False
    printA = False
    saveA = False

angle = 'mili-rad'
length = 'nm'
volt = 'mili-volt'

if angle == 'mili-rad':
    X_diff = linspace(0,1,N_plot)

if length == 'micron':
    a = 85 # micron  
    d = .406E2 # micron
    hc = 1.2398 # eV micron
    red,e_red = .670,.005 # micon
    green,e_green = .546,.005 # micron
if length == 'nm':
    a = 85000 # nm
    d = 406 # nm
    hc = 1239.8 # eV nm
    red,e_red = 670,5 # nm
    green,e_green = 546,5 # nm



# Error Tracking
e_cnt = .5 # num
dt,e_dt = .3,.05 # s
red,e_red = 670,5 # nm
green,e_green = 546,5 # nm
L,e_L = 4.06E-4, .01 # meter
e_A,e_V = .005,.005 # radian,volt
noise_cnt,e_noise_cnt = 39.80,6.421 # num
noise_V,e_noise_V = .0221/1000,.003564/1000 # V
V0,e_V0 = .2272/1000,.01370/1000 # V
Cmax,e_Cmax = 369.5,24.67 # num
    

# Program Functions
def OpenCSV(TITLE,argN,argA):
    DATA = []
    with open(TITLE, newline='', encoding='utf-8') as d:
        reader = csv.reader(d)
        for row in reader:
            DATA.append(row)
    if argA != 0:
        return array(DATA[argN:]).astype(float)
    else:
        return DATA[argN:]

#############################################################3
    
def WriteCSV(TITLE,ARRAY):
    with open(TITLE, 'w', newline='') as fp:
        write = csv.writer(fp, delimiter=',')
        write.writerows(ARRAY)
    return

##########################################################

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()   
    return array[idx]

##############################################################

def Stats(ARRAY):
    ARRAY = array(ARRAY).astype(float)
    STD = std(ARRAY)
    VAR = var(ARRAY)
    MEAN = average(ARRAY)
    return [STD,VAR,MEAN]

################################################################################
    
def BinW(ARRAY,N):
    ARRAY = array(ARRAY).astype(float)
    MIN = np.min(ARRAY)
    MAX = np.max(ARRAY)
    return (MAX-MIN)/N

############################################################################

def n_std(z):
    return erf(z/sqrt(2))    

################################################################################

def BOX(axis,X,Y,color1,color2,kind):
    X,Y = array(X),array(Y)
    Yav = []
    
    Y2 = np.sort(Y)
    yav = average(Y2)
    
    for x in X:
        Yav.append(yav)
        
    xmin = np.min(X)
    width = np.max(X)-np.min(X)
    
    if abs(max(Y2)) > abs(min(Y2)):
        ymax = max(Y2)
    else:
        ymax = abs(min(Y2))
    
    def PercentLimits():
        dy = .1
        hw = dy
        #target = goal
        PER,HW = [],[]
        
        while hw < ymax:
            LL,UL = yav-hw,yav+hw
            Li,Ui = (np.abs(Y2-LL)).argmin(), (np.abs(Y2-UL)).argmin()
            per = len(Y2[Li:Ui])/len(Y2)
            PER.append(per)
            HW.append(hw)
            hw += dy
        
        PER,HW = array(PER),array(HW)
        return PER,HW
    PER,HW = PercentLimits()
    
    if kind == 'sig1':
        hw = std(Y)
        ymin = yav-hw
        height = 2 * hw
        
        #ind = (np.abs(HW-hw)).argmin()
        per = .65 #PER[ind]
        
        #LL,UL = yav-hw,yav+hw
        #Li,Ui = (np.abs(Y2-LL)).argmin(), (np.abs(Y2-UL)).argmin()
        #per2 = len(Y2[Li:Ui])/len(Y2)
        delta_per = .001 #(per-per2)
    
    if kind == 'sig2':
        per = n_std(2)
        ind = (np.abs(PER-per)).argmin()
        hw = round(HW[ind],2)
        
        ymin = yav-hw
        height = 2 * hw
        
        LL,UL = yav-hw,yav+hw
        Li,Ui = (np.abs(Y2-LL)).argmin(), (np.abs(Y2-UL)).argmin()
        per2 = len(Y2[Li:Ui])/len(Y2)
        delta_per = per-per2
    
    if kind == '68':
        per = .68
        ind = (np.abs(PER-per)).argmin()
        hw = HW[ind]
        
        ymin = yav-hw
        height = 2 * hw
        
        LL,UL = yav-hw,yav+hw
        Li,Ui = (np.abs(Y2-LL)).argmin(), (np.abs(Y2-UL)).argmin()
        per2 = len(Y2[Li:Ui])/len(Y2)
        delta_per = per-per2
            
    if kind == '95':
        per = .95
        ind = (np.abs(PER-per)).argmin()
        hw = HW[ind]
        
        ymin = yav-hw
        height = 2 * hw
        
        LL,UL = yav-hw,yav+hw
        Li,Ui = (np.abs(Y2-LL)).argmin(), (np.abs(Y2-UL)).argmin()
        per2 = len(Y2[Li:Ui])/len(Y2)
        delta_per = per-per2
    
    print("%s statistics: half-width = %s , per = %s, delta_per = %s" % (kind,hw,per,delta_per))
        
    p = patches.Rectangle((xmin,ymin),width,height,facecolor=color1,alpha=.6)
    plot1 = axis.add_patch(p)
    plot2 = axis.plot(X,Yav,color2)
    return plot1,plot2,per,hw

##################################################################################

def Leg(param,COLOR,LABEL):
    if param == 'data':
        return mlines.Line2D([], [], color=COLOR, marker='s', markersize=10, label=LABEL,linewidth=5)
    else:
        return mlines.Line2D([], [], color=COLOR, marker="", markersize=10, label=LABEL,linewidth=5)
def Legend(AXIS,LegArray,LOC):
    return AXIS.legend(handles = LegArray, loc = LOC, numpoints = 1,fontsize=20) 

##################################################################

def XY(typeA,TABLE,YCorrection): # average array
    X,Y = [],[]
    if typeA == 'diff':
        for row in TABLE:
            X.append(row[1]) # mili-rad
            Y.append((row[2]+YCorrection)) # mili-volt
    else:
        for row in TABLE:
            X.append(row[0]) # mili-rad
            Y.append(row[1]+YCorrection) # mili-volt             
    return array(X).astype(float),array(Y).astype(float)
    
########################################################################

def AxisPlot(axisA,LOC,X,Y,COLOR,XLABEL,YLABEL,HANDLES,plotA):
    axisA = plt.subplot(LOC[0],LOC[1],LOC[2])
    axisA.set_xlabel(XLABEL,fontsize=labS)
    axisA.set_ylabel(YLABEL,fontsize=labS)
    Legend(axisA,HANDLES,1)
    
    if plotA == 'plot':     
        axisA.plot(X,Y,COLOR,linewidth=.5)
    return axisA

########################################################################

def CntToVolt(C,dt,waveL,VCorrection): # ( # , s , nm )
    C = array(C).astype(float)
    def V(c):
        return (c * hc * dt)/(waveL * 1000) # mV
    return  V(C) + VCorrection

###########################################################################

def deltaV(C,ec,dt,edt,wL,ewL):
    Y = []
    def dV(c):
        try:
            result = (hc/1000) * sqrt( (dt*ec/wL)**2 + (c*edt/wL)**2 + (c*dt/(ewL**2))**2 )
        except TypeError:
            result = -(hc/1000) * sqrt( abs((dt*ec/wL)**2 + (c*edt/wL)**2 + (c*dt/(ewL**2))**2 ) )
        return result
    for c in C:
        Y.append(dV(c))
    Y = array(Y).astype(float)
    return Y
        
 ####################################################################   

def SS(x,par):
    A,B,C,D = par[0],par[1],par[2],par[3]
    theta = (x - B)
    alpha = (pi*C/D) * (theta)
    return A * (sin(alpha))/alpha
    
def DS(x,par):
    A,B,C,D,E = par[0],par[1],par[2],par[3],par[4]
    theta = (x-B)
    alpha = (pi*C/D) * theta
    beta = (pi*E/C) * theta
    return A * cos(beta)**2 * (sin(alpha)/alpha)**2

#################################################################################

def SSFit1(slit,X,Y,par):
    ymax,shift,a,b = par[0],par[1],par[2],par[3] 

    def fitfunc(p,x):
        theta = (x-shift)
        return p[0] * sin((pi*p[1]*theta/p[2])) / ((pi*p[1]*theta/p[2]))
    def errfunc(p, x, y):
        return y - fitfunc(p, x)

    xdata=array(X)
    ydata=array(Y)
    
    if slit == 'A':
        cycles = 5000
    if slit == 'B':
        cycles = 8000
    
    qout,success = optimize.leastsq(errfunc, [ymax,a,b] , args=(xdata, ydata),maxfev=cycles)

    out = qout[:]
    out[0] = qout[0]
    out[1] = qout[1]
    out[2] = qout[2]
    fitParams = [out[0],shift,out[1],out[2]]
    #if success == True:
    #    print("success")
    #else:
    #    print("OMFG! FAILED!!!")
    
    def FUNC(x):
        if slit == 'A':
            theta = (x - shift)
        if slit == 'B':
            theta = (x - shift)
        alpha = (pi*out[1]/out[2]) * theta
        return out[0] * (sin(alpha))/alpha
    
    Y_out = FUNC(X_diff)
    return Y_out,fitParams

################################################################################

def SSFit2(slit,X,Y,par):
    ymax,shift,a,b = par[0],par[1],par[2],par[3]    

    def fitfunc(p,x):
        if slit == 'A':
            theta = (x - shift)
        if slit == 'B':
            theta = (x - shift)
        return ymax * sin( (pi*p[0]*theta/p[1]) ) / ( (pi*p[0]*theta/p[1]) )
    def errfunc(p, x, y):
        return y - fitfunc(p, x)

    xdata=array(X)
    ydata=array(Y)
    
    if slit == 'A':
        cycles = 5000
    if slit == 'B':
        cycles = 5000
    
    qout,success = optimize.leastsq(errfunc, [a,b] , args=(xdata, ydata),maxfev=cycles)

    out = qout[:]
    out[0] = qout[0]
    out[1] = qout[1]
    fitParams = [ymax,shift,out[0],out[1]]
    #if success == True:
    #    print("success")
    #else:
    #    print("OMFG! FAILED!!!")
    
    def FUNC(x):
        if slit == 'A':
            theta = (x - shift)
        if slit == 'B':
            theta = (x - shift)
        alpha = (pi*out[0]/out[1]) * theta
        return ymax * (sin(alpha))/alpha
    
    Y_out = FUNC(X_diff)
    return Y_out,fitParams
    
#######################################################################################    
    
def DSFit(X,Y,par):
    A,B,C,D,E = par[0],par[1],par[2],par[3],par[4]   
    
    def fitfunc(p,x):
        theta = (x-B)
        alpha = (pi*p[1]*theta)/p[2]
        beta = (pi*p[3]*theta)/p[2]
        return 4*p[0] * cos(beta)**2 * (sin(alpha)/alpha)**2
    def errfunc(p, x, y):
        return y - fitfunc(p, x)

    xdata=array(X)
    ydata=array(Y)
    
    qout,success = optimize.leastsq(errfunc, [A,C,D,E] , args=(xdata, ydata),maxfev=6000)

    out = qout[:]
    out[0] = qout[0] # V0
    out[1] = qout[1] # a
    out[2] = qout[2] # waveL
    out[3] = qout[3] # d
    fitParams = [out[0],B,out[1],out[2],out[3]]
    #print("")
    #print("fit params:")
    #print(fitParams)
    #print("")
      
    def FUNC(x):
        theta = (x-B)
        alpha = (pi*out[1]/out[2]) * theta
        beta = (pi*out[3]/out[2]) * theta
        return 4*out[0] * cos(beta)**2 * (sin(alpha)/alpha)**2
    
    Y_out = FUNC(X_diff)
    return Y_out,fitParams

####################################################################################

def PrintParams(slit,ARRAY): 
    ARRAY[0] = round(ARRAY[0],2)
    R1 = (ARRAY[2]/ARRAY[3])*100
        
    if len(ARRAY) == 4:
        a = round(R1*red,3)
        print("a/wavelength = %s" % round(R1,3))
        print("slit width 'a' = R1*lambda = %s nm" % a)
        #print("%s fitting params1 : V_0 = %s %s, shift = %s %s, a = %s %s, lambda = %s %s" % (slit,ARRAY[0],volt,ARRAY[1],angle,ARRAY[2],length,ARRAY[3],length))
    
    if len(ARRAY) == 5:
        R2 = (ARRAY[4]/ARRAY[2]) * 100
        a = round(R1*red,3)
        d = round(R2*red,3)
        print("a/wavelength = %s" % a)
        print("slit width 'a' = R1 * lambda = %s nm" % a)
        print("d = %s nm" % d)
        #print("%s fitting params1 : V_0 = %s %s, shift = %s %s, a = %s %s, lambda = %s %s, d = %s %s" % (slit,ARRAY[0],volt,ARRAY[1],angle,ARRAY[2],length,ARRAY[3],length,ARRAY[4],length))
    return
    
def findA(a1,a2,a3):
    As = np.array([a1,a2,a3]).astype(float)
    print(As)
    return average(As),std(As)

#####################################################################################

def Hist_Plot(figN,X,Y,COLOR,BINWIDTH):
    fig = plt.figure(figN,figsize=(30,15))
    nullfmt = NullFormatter()
    
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    
    #axScatter.legend(handles = HANDLES, loc = 4,numpoints = 1)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # the scatter plot:
    axScatter.scatter(X,Y,color = 'k')
    BOX(axScatter,X,Y,COLOR,'r','sig1')
    axScatter.set_yticks([average(Y),average(Y)-std(Y),average(Y)+std(Y)])

    # now determine nice limits by hand:
    binwidth = BINWIDTH
    xymax = np.max([np.max(np.fabs(X)), np.max(np.fabs(Y))])
    lim = (int(xymax/binwidth) + 1) * binwidth

    axScatter.set_xlim((min(X), max(X)))
    axScatter.set_ylim((max(Y),min(Y)))
    
    bins = np.arange(-lim, lim + binwidth, binwidth)
    axHistx.hist(X, bins=bins, color=COLOR)
    axHisty.hist(Y, bins=bins, orientation='horizontal', color=COLOR)
    
    axHistx.set_xlim(axScatter.get_xlim())
    axHisty.set_ylim(axScatter.get_ylim())
    
    axScatter.set_xlabel('t (s)',fontsize=fontS)
    axScatter.set_ylabel('Counts',fontsize=fontS)
    
    return fig

######################################################################################
  
def normalize(v):
    norm=np.linalg.norm(v)
    if norm==0: 
       return v
    return v/norm
    
################################################################################

def PLOTS(figN,binsize,TITLE,fitA,par):
    TABLE = OpenCSV(TITLE,1,1)
    
    Y,L = [],[]
    for row in TABLE:
        Y.append(row[1])
        
    num_bins = ((max(Y)-min(Y))/binsize)
    
    cnt_cont = linspace(min(Y),max(Y),N_plot)
    cnt_bin = linspace(min(Y),max(Y),num_bins)
    
    COUNT = 0
    for c in cnt_bin[:num_bins-1]:
        i = (np.abs(cnt_bin-c)).argmin()
        count = 0
        length = 0
        for row in TABLE:
            y = row[1]
            if y >= cnt_bin[i] and y < cnt_bin[i+1]:
            #if abs(c-y) <= binsize/2:
                count += y
                length += 1
                COUNT += 1
        if length != 0:
            L.append(count/length)
        if length == 0:
            L.append(0)
            
    L = array(L)#.astype(int)
    L_N = normalize(L)
    cnt_av = np.average(Y)
    Lambda = cnt_av

    print("COUNT = %s" % COUNT)
    print("samples dropped = %s" % (500-COUNT))
    print("Lambda = %s" % Lambda)
    print("Normed: k_var = %s , k_av = %s" % (np.var(L_N),np.average(L_N)))
    print("Not Normed: k_var = %s , k_av = %s" % (np.var(L),np.average(L)))
    
    def histogram(ARRAY):
        
        sigma = np.std(ARRAY)  # standard deviation of distribution
        x = sigma * Y
        
        n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
        plt.subplots_adjust(left=0.15)
        return
   
    def poisson(C,shift):
        def P(c):
            x = c
            return Lambda**(x) * np.exp(-Lambda) / factorial(x)
        
        C = C-shift
        
        return P(C)
    YP = poisson(cnt_cont,0)
    
    def gaussian(C,F): # count and frequency
        
        def FIT(X,Y):
            def fitfunc(p,x):
                one = 1/np.sqrt(2*p[0]*pi)
                two_1 = -(x-p[1])**2
                two_2 = 2*p[0]
                two = np.exp(two_1/two_2)
                return one*two
            def errfunc(p,x,y):
                y-fitfunc(p,x)
        
            xdata=array(X)
            ydata=array(Y)

            qout,success = optimize.leastsq(errfunc, [par[0],par[1]], args=(xdata, ydata),maxfev=5000)

            out = qout[:]
            out[0] = qout[0]
            out[1] = qout[1]
            return out[0],out[1]
       
        if fitA == True:
            var,av = FIT(cnt_bin[:num_bins-1],L_N)
        else:
            var,av = par[0],par[1]
        
            
        def FUNC(x,variance,mean):
            one = 1/sqrt(2*variance*pi)
            two_1 = -(x-mean)**2
            two_2 = 2*variance
            two = np.exp(two_1/two_2)
            return one*two
        return FUNC(cnt_cont,var,av)
    YG = gaussian(cnt_cont,L_N)
    
    HL = Leg('data','g','Data')
    PL = Leg(0,'b','Poisson')
    GL = Leg(0,'r','Gaussian')
    
    fig = plt.figure(figN,figsize=figS)
    plt.xlim([min(Y),max(Y)])
    #plt.xticks(cnt_bin)
    plt.xlabel('Counts',fontsize=fontS)
    plt.ylabel('Count Frequency',fontsize=fontS)
    Legend(plt,[HL,GL,PL],1)
    plt.grid(False)
    
    #plt.hist(K,cnt_bin,facecolor='g',normed=True)
    histogram(Y)
    plt.plot(cnt_cont,YP,'b-',linewidth=2)
    #plt.scatter(cnt_bin[:num_bins-1],poisson(K_N),color='b',s=30)
    plt.plot(cnt_cont,YG,'r-',linewidth=2)
    
    return fig

    

  

###################################################################################


####################################################################################


#######################################################################################



 


