import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import optimize
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import csv

###############################################################################

# CONSTANTS

path = 'C:\\Users\\Jacob\\Google Drive\\School\\Physics Degree\\Current Classes\\PHY 334 Advanced Lat 1\\Compton\\'
CSVpath = path + 'CSV\\'
NPpath = path + 'NP\\'
Gpath = path + 'Graphs\\'
size = (30,15)
LW = .5
FS = 20

binmin,binmax = 56,98
Ba,Cs = 356,661 # keV
A,B = -50.666666666666686, 7.261904761904762
E0e = 511 # rest energy of electron keV

N_A = 6.022E23 # 1/mols
#conversions
mCi = 3.7e7 # decays/s
inch = 2.54 # cm
gram = 0.03706237946850917 # mols


###############################################################################




###############################################################################

# CLASSES

class parameters:
    def __init__ (self,LL,xmin,xmax,ymin,ymax,xlabel,ylabel,convert):
        self.LL = LL
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.convert = convert
                   #(LL,xmin, xmax, ymin, ymax, xlabel,             ylabel,             convert)
params0 = parameters(1, 20,   110,  0,    1600, 'Channel #',        'CPS',              False) # bins
params1 = parameters(1, 100,  800,  0,    1600, 'keV',              'CPS',              True) # keV
params2 = parameters(1, 0,    150,  0,    4000, 'Channel #',        'CPS',              False) # bins
params3 = parameters(1, 100,  800,  0,    4000, 'keV',              'CPS',              True) # keV
params4 = parameters(1, 100,  800,  0,    600,  'keV',              'CPS',              True) # keV
params5 = parameters(1, 0,    1,    0,    .005, '1-cos($\\theta$)', '(1/E) [keV$^-1$]', True)

class angle:
    def __init__ (self,index,deg,CSV,axis,label):
        self.index = index
        self.deg = deg
        self.CSV = CSVpath + CSV + '.csv'
        self.axis = axis
        self.label = label

calibration = angle(0,0,'Calibration',111,'Calibration')
        
        #   (i,deg, CSV,    axis,  label)
in20 = angle(0,20,  'in20', 241,   '20$^o$: w/')
in30 = angle(1,30,  'in30', 242,   '30$^o$: w/')
in40 = angle(2,40,  'in40', 243,   '40$^o$: w/')
in50 = angle(3,50,  'in50', 244,   '50$^o$: w/')
in60 = angle(4,60,  'in60', 245,   '60$^o$: w/')
in70 = angle(5,70,  'in70', 246,   '70$^o$: w/')
in80 = angle(6,80,  'in80', 247,   '80$^o$: w/')
in90 = angle(7,90,  'in90', 248,   '90$^o$: w/')
IN = [in20,in30,in40,in50,in60,in70,in80,in90]

        #    (i,deg, CSV,    axis,   label)
out20 = angle(0,20,  'out20', 241,   '20$^o$: w/o')
out30 = angle(1,30,  'out30', 242,   '30$^o$: w/o')
out40 = angle(2,40,  'out40', 243,   '40$^o$: w/o')
out50 = angle(3,50,  'out50', 244,   '50$^o$: w/o')
out60 = angle(4,60,  'out60', 245,   '60$^o$: w/o')
out70 = angle(5,70,  'out70', 246,   '70$^o$: w/o')
out80 = angle(6,80,  'out80', 247,   '80$^o$: w/o')
out90 = angle(7,90,  'out90', 248,   '90$^o$: w/o')
OUT = [out20,out30,out40,out50,out60,out70,out80,out90]

INOUT = IN + OUT

class diff:
    def __init__ (self,index,deg,peakx,peaky,NP,QE,PtT,width,fit_manual,fit_alg,Gwidth,CS):
        self.index = index
        self.deg = deg
        self.peakx = peakx
        self.peaky = peaky
        self.NP = NPpath + NP
        self.inA = IN[self.index]
        self.outA = OUT[self.index]
        self.deltaE = Cs-self.peakx
        self.QE = QE
        self.PtT = PtT
        self.width = width
        self.ROI_min = peakx - width
        self.ROI_max = peakx + width
        self.fit_manual = fit_manual # [ height , average , sigma]
        self.fit_alg = fit_alg
        self.Gwidth = Gwidth
        self.LLg = self.fit_alg[1]-self.Gwidth
        self.ULg = self.fit_alg[1]+self.Gwidth
        self.CS = CS # Cross section 1E-27 cm^2
         #(i, deg, xpeak, ypeak, NP,    QE,   PtT, width, fit_manual,   fit_alg,                                      width,        CS)
d20 = diff(0, 20,  602,   440,   'd20', .865, .47, 100,   [585,602,30], [527.52806165,  593.2137932,    39.03494422], 143.2137932,  1)
d30 = diff(1, 30,  552,   329,   'd30', .890, .50, 100,   [450,552,30], [394.72638306,  546.83202491,   43.49193105], 154.83202491, 1)
d40 = diff(2, 40,  479,   269,   'd40', .930, .53, 80,    [362,479,30], [304.04510625,  492.33180262,   46.31160992], 166.33180262, 1)
d50 = diff(3, 50,  428,   221,   'd50', .950, .55, 90,    [279,428,30], [247.24221962,  428.93662367,   39.6965504 ], 145.93662367, 1)
d60 = diff(4, 60,  392,   198,   'd60', .960, .57, 80,    [294,392,30], [243.50208966,  388.3122186,    40.02825416], 141.3122186,  1)
d70 = diff(5, 70,  341,   204,   'd70', .980, .61, 90,    [298,341,30], [225.48489262,  335.42059097,   55.13930632], 197.42059097, 1)
d80 = diff(6, 80,  290,   245,   'd80', .990, .65, 60,    [334,290,30], [273.21748153,  299.01899898,   31.19974875], 111.01899898, 1)
d90 = diff(7, 90,  268,   208,   'd90', .995, .69, 60,    [295,268,30], [214.71576835,  273.9431573,    35.34852891], 128.9431573,  1)
DIFF = [d20,d30,d40,d50,d60,d70,d80,d90]

###############################################################################




###############################################################################

# MISC FUNCTIONS

def OpenCSV(TITLE,LL,arrayA): # Imports Table from csv file in python file folder
    # LL,UL give the interval of the table rows to import
    DATA = []
    with open(TITLE, newline='', encoding='utf-8') as d:
        reader = csv.reader(d)
        for row in reader:
            DATA.append(row)
    if arrayA != 0:
        return np.array(DATA[LL:]).astype(arrayA)
    else:
        return DATA[LL:]

def AB(x1,x1p,x2,x2p):
    B = (x2p-x1p)/(x2-x1)
    A = x1p - B*x1
    print(A,B)
    return A,B
#AB(binmin,Ba,binmax,Cs)
    
def BinToKeV(A,B,x): # converts "channel" to "keV"
    x_prime = A + B*x
    return x_prime # keV

def KeVToBin(A,B,x_keV): # converts "keV" to channel
    x_bin = (1/B) * (x_keV - A)
    if x_bin < 0:
        return 0
    else:
        return int(x_bin)

def XY(data,convertA):
    X,Y = [],[]
    for row in data:
        if convertA == True:
            X.append(BinToKeV(A,B,row[0])) # keV
        else:
            X.append(row[0])
        Y.append(row[1]) # CPS
    X = np.array(X).astype(int)
    Y = np.array(Y).astype(int)
    return X,Y

def n_std(z):
    return sp.special.erf(z/sp.sqrt(2))

def PLOT(case,params,data_list,loc):
    fig = plt.figure(case,figsize=size)
    plt.xlabel(params.xlabel,size=FS)
    plt.ylabel(params.ylabel,size=FS)
    plt.xlim([params.xmin,params.xmax])
    plt.ylim([params.ymin,params.ymax])
    
    HANDLES = []
    for d in data_list:
        data = OpenCSV(d.CSV,params.LL,int)
        X,Y = XY(data,params.convert)
        l = plt.plot(X,Y, label=d.label)
        HANDLES.append(l)
    
    if len(data_list)>1:
        plt.legend(loc=loc,fontsize=20)
        
    return fig

def PLOT_axis(case,params,data_list):
    fig = plt.figure(case,figsize=size)
    
    AX,HANDLES = [],[]
    for d in data_list:
        ax = fig.add_subplot(d.axis)
        AX.append(ax)
        ax.set_xlabel(params.xlabel,size=FS)
        ax.set_ylabel(params.ylabel,size=FS)
        ax.set_xlim([params.xmin,params.xmax])
        ax.set_ylim([params.ymin,params.ymax])
        data = OpenCSV(d.CSV,params.LL,int)
        X,Y = XY(data,params.convert)
        l = plt.plot(X,Y, label=d.label)
        HANDLES.append(l)
    
    for a in AX:
        a.legend(loc=2,fontsize=20)
    
    plt.tight_layout()
    
    return fig,AX
        



#def plot_axis()
    
    
    

###############################################################################




###############################################################################

# LAB FUNCTIONS

def Eq2(x,E_0,E): # x = 1-cos(theta) , E_0 = mc^2 of electron
    return (x/E_0) + (1/E)

def diff_crossSection(yeild,ret='one'):
    # I 
    Ro = 5.46 * mCi # decays/s
    yr = 21.5 # yrs after calibration
    R = Ro * np.exp(-yr/43.28) # decay rate
    #print(R/1e8)
    
    # II
    r_ts = 6 * inch # distance from source to target in cm
    A = 4 * np.pi * (r_ts**2) # Area of sphere at target distance from source in cm^2
    #print(A)
    
    # III
    I_inc = R/A # incident intensity at target in (#/s)/cm^2
    #print(I_inc)
    
    # IV
    d_t = 1.4 # diameter of target in cm
    r_t = d_t/2 # radius of target in cm
    h_t = 1.9 # height of target in cm
    V_t = (np.pi*r_t**2)*h_t # volume of target in cm^3
    rho_al = 2.7 # g/cm^3
    N_atoms = V_t * rho_al * gram * N_A # number of atoms
    n_electrons = 13 # electrons per Al atom
    N_electrons = N_atoms * n_electrons # number of elctrons
    #print(N_electrons)
    
    # V
    d_d = 1.75 * inch # diameter of detector in cm
    A_d = np.pi * d_d**2 * (1/4)
    r_dt = 14 * inch # distance from target to detector in cm
    d_Omega = A_d/r_dt**2 # differential solid angle in sr
    #print(d_Omega)
    #print(1/(d_Omega*N_electrons*I_inc*150))
    
    e1 = 7.262/(d_Omega*N_electrons*I_inc)
    e2 = ((np.pi*d_d*.127)/(2*r_dt**2) + (np.pi*d_d**2*.127)/(2*r_dt**3)) * (yeild)/(N_electrons*I_inc*d_Omega**2)
    e3 = ((2*N_electrons*.127)/r_t + (N_electrons*.127)/h_t) * (yeild)/(d_Omega*I_inc*N_electrons**2)
    unc = np.sqrt(e1**2 + e2**2 + e3**2 + 7.262**2)
    
    # VI
    if ret == 'one':
        return (yeild)/(d_Omega * N_electrons * I_inc*150) # cm^2
    else:
        return unc
    
def KleinNishina(theta):
    theta = theta * (np.pi/180) # radians
    r0 = 2.82e-13 # cm
    gamma = Cs/E0e
    a1 = 1 - np.cos(theta)
    a2 = 1 + np.cos(theta)**2
    one = r0**2
    two = a2/2
    three = (1 + gamma*a1)**(-2)
    four_a = (gamma * a1)**2
    four_b = a2 * (1 + gamma*a1)
    four = 1 + (four_a/four_b)
    ans = one * two * three * four
    return ans*1e27 # cm^2

def gaus(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))

def PeakTotal(X,Y,Y_g): # insert total and gaussian arrays
    sum1,sum2 = 0,0
    for i in range(len(X)):
        sum1 += abs(Y[i])
        sum2 += abs(Y_g[i])
    return sum2/sum1    
    
###############################################################################




###############################################################################

# FITTING FUNCTIONS

def GaussFit(X,Y,d):
    
    LL = (np.abs(X-d.ROI_min)).argmin()
    UL = (np.abs(X-d.ROI_max)).argmin() + 1
    X2 = ar(X[LL:UL])
    Y2 = ar(Y[LL:UL])

    def gaus(x,a,x0,sigma):
        return a*exp(-(x-x0)**2/(2*sigma**2))

    popt,pcov = curve_fit(gaus,X2,Y2,p0=d.fit_manual)
    return LL,UL,popt

def ComptonFit(X,Y):
    
    def fitfunc(p,x):
        return (x/p[0]) + (1/p[1])
    def errfunc(p, x, y):
        return y - fitfunc(p, x)

    xdata=np.array(X)
    ydata=np.array(Y)
    
    qout,success = optimize.leastsq(errfunc, [511,Cs] , args=(xdata, ydata),maxfev=5000)

    out = qout[:]
    out[0] = qout[0]
    out[1] = qout[1]
    
    X2 = np.linspace(0,1,1000)
    Y2 = Eq2(X2,out[0],out[1])
    return X2,Y2,out[0],out[1]
    
    
    
    
    
    
    
    
    
    