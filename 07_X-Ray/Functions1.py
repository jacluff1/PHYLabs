import numpy as np
import matplotlib.lines as mlines
import csv

############################################################

# CONSTANTS

############################################################

R = 13.6 # ev, boar hydrogen ionizizing energy
R_F = 3.29E15 # Rydberg fundamental frequecy Hz
h = 4.135667662E-18 # Planck's constant keV s

C_alpha = (R_F) * (3/4) * h # expected value for eq 2
sigma_alpha = 1

size = (30,14) # figsize
FS = 20 # font size for labels
lw = 1
ms = 20
Xlabel = 'keV'
Ylabel = 'Counts'

#A,B = 0.372,0.0283 # old
A,B = -.1330,.03027 # new
dA,dB = .00006,.000001

Mo_I0 = 1429
rho_Zr = 6.44 #g cm^3

############################################################

# CLASSES

############################################################

class Carousel:
    def __init__ (self,index,Z,name,title,C,sigma,err_mosley,params,ROI,loc,alpha,beta,  SO_ROI,doublet,err_doub,doub_axis,D_ROI,doubFit):
        self.index = index # index in the samples list
        self.Z = Z # atomic number
        self.name = name # element name
        self.title = title # csv title
        self.C = C # mosley fitting parameter 1
        self.sigma = sigma # mosley fitting parameter 2
        self.err_mosley = err_mosley # error from mosley fit
        self.params = params # [loc,color,style]
        self.ROI = ROI # keV
        self.loc = loc # legend location
        self.alpha = alpha # K_alpha keV
        self.beta = beta # K_beta keV
        self.SO_ROI = SO_ROI
        self.doublet = doublet
        self.err_doub = err_doub
        self.doub_axis = doub_axis
        self.D_ROI = D_ROI
        self.doubFit = doubFit
#           index,Z,name,title,  C,                sigma,           err_mosley,       params,          ROI,loc,alpha,beta,     SO_ROI,   doublet,      err_doub, doub_axis,D_ROI,     doubFit)
Cu = Carousel(0,29,'Cu','Cu.csv',0.0102272329008 , 0.944444444444 , 2.17700431904e-07,[231,'r','rs',1],[5,10],2,8.05,8.9,      [8.5,12], 'no',         'no',     'no', 'no',          'no')
Rb = Carousel(1,37,'Rb','Rb.csv',0.0103579845799 , 1.04554554555 , 2.77797467518e-06,[232,'g','gs',1],[10,20],1,13.39,14.96,   [14,18] , 'no',         'no',     'no', 'no' ,         'no')
Mo = Carousel(2,42,'Mo','Mo.csv',0.0103130386903 , 0.83033033033 , 6.6285794631e-07,[233,'c','cs',1],[10,25],2,17.443,19.651,  [19,23] , 'no',         'no',     'no', 'no',          'no')
Ag = Carousel(3,47,'Ag','Ag.csv',0.0104111024496 , 0.864364364364 , 9.65495324579e-07,[234,'b','cs',1],[15,30],2,22.16,24.94,  [24,27] , [25.45,120] , [.005,5], 311,  [25.35,25.47], [4.011, 0.1765])
Ba = Carousel(4,56,'Ba','Ba.csv',0.0106971217477 , 1.14364364364 , 9.08419715984e-06,[235,'m','ms',1],[25,40],2,32.19,36.38,   [35,39] , [37.25,15] , [.005,.5], 312,  [37.1,37.4],   [3.969, 0.2335])
Tb = Carousel(5,65,'Tb','Tb.csv',0.0109095932262 , 1.15465465465 , 9.28880223938e-06,[236,'k','ks',1],[35,60],1,44.47,50.39,   [49,54] , [51.7,260] , [.05,5],   313,  [51.6,51.75],  [3.925, 0.3958])
Samples = [Cu,Rb,Mo,Ag,Ba,Tb]

class MO_thickness:
    def __init__ (self,index,name,csv, thick , lam , mu_E , plot_params):
        self.index = index
        self.name = name
        self.csv = csv
        self.thick = thick # in
        self.lam = lam # [ lambda , attenuation , I_attenuated , error_per ]
        self.mu_E = mu_E # g^-1 cm^-4
        self.plot_params = plot_params # [ subplot , legend , color ]
T1 = MO_thickness(0,'.01 inch','001.csv' , .01 , [0.008108, 416.3, 1013, 14.3, 1.412] , [7.540, 33780] , [231,2,'r'])
T2 = MO_thickness(1,'.02 inch','002.csv' , .02 , [0.02223, 581.0, 848.0, .002525, 0.0002978] , [2.750, .7934] , [232,2,'orange'])
T3 = MO_thickness(2,'.03 inch','003.csv' , .03 , [0.05608, 837.0, 592.0, .03679, 0.006214] , [1.090, 1.816] , [233,2,'g'])
T4 = MO_thickness(3,'04 inch','004.csv' , .04 , [0.08995, 916.0, 513.0, .001361, 0.0002653] , [0.6796, .02612] , [234,2,'c'])
T5 = MO_thickness(4,'.05 inch','005.csv' , .05 , [0.1509, 1026, 403.0, .03250, 0.008066] , [0.4051, .2216] , [235,2,'b'])
T6 = MO_thickness(5,'.06 inch','006.csv' , .06 , [0.09231, 746.0, 683.0, .008950, 0.001310] , [0.6622, .1631] , [236,2,'m'])
thickness = [T1,T2,T3,T4,T5,T6]

class WheelT:
    def __init__ (self,index,name,csv,plot,E_0,E_A):
        self.index = index
        self.name = name
        self.csv = csv
        self.plot = plot # [ subplot , legend , color ]
        self.E_0 = E_0
        self.E_A = E_A
        #    i,'name',CSV,   , legend,    
CuT = WheelT(0,'Cu','Cu.csv' , [231,2,'r'] , 1 , 1)
RbT = WheelT(1,'Rb','Rb.csv' , [232,2,'g'] , 1 , 1)
MoT = WheelT(2,'Mo','Mo.csv' , [233,2,'c'] , 1 , 1)
AgT = WheelT(2,'Ag','Ag.csv' , [234,2,'b'] , 1 , 1)
BaT = WheelT(3,'Ba','Ba.csv' , [235,2,'m'] , 1 , 1)
TbT = WheelT(4,'Tb','Tb.csv' , [236,2,'k'] , 1 , 1)
Wheel = [CuT,RbT,MoT,AgT,BaT,TbT]



#############################################################

# MISC FUNCTIONS

############################################################

def OpenCSV(TITLE,LL,UL,numpyA): # Imports Table from csv file in python file folder
    # LL,UL give the interval of the table rows to import
    # numpyA gives data type: "float" , "int" , or "str"
    DATA = []
    with open(TITLE, newline='', encoding='utf-8') as d:
        reader = csv.reader(d)
        for row in reader:
            DATA.append(row)
    if numpyA == float or numpyA == int:
        return np.array(DATA[LL:UL]).astype(numpyA)
    else:
        return DATA[LL:UL]

def WriteCSV(TITLE,ARRAY):
    with open(TITLE, 'w', newline='') as fp:
        write = csv.writer(fp, delimiter=',')
        write.writerows(ARRAY)
    return

def BinToKeV(A,B,x): # converts "channel" to "keV"
    x_prime = A + B*x
    return x_prime # keV

def KeVToBin(A,B,x_keV): # converts "keV" to channel
    x_bin = (1/B) * (x_keV - A)
    if x_bin < 0:
        return 0
    else:
        return int(x_bin)

def XY(convertA,TABLE,i_channel,i_count): # Converts table data into X and Y arrays
    X,Y = [],[]
    for row in TABLE:
        if convertA == False:
            X.append(row[i_channel])
        else:
            keV = round( BinToKeV(A,B,row[i_channel]),2) # converts "channel" to "keV"
            X.append(keV)
        Y.append(row[i_count])
    X,Y = np.array(X).astype(float),np.array(Y).astype(int)
    return X,Y

def Leg(param,COLOR,LABEL):
    if param == 'data':
        return mlines.Line2D([], [], color=COLOR, marker='s', markersize=10, label=LABEL,linewidth=5)
    else:
        return mlines.Line2D([], [], color=COLOR, marker="", markersize=10, label=LABEL,linewidth=5)

def Legend(AXIS,LegArray,LOC):
    return AXIS.legend(handles = LegArray, loc = LOC, numpoints = 1,fontsize=20) 

def plotAxis(fig,axis,plotA,X,Y,el,xticks):
        
    axis = fig.add_subplot(el.params[0])
    axis.set_xlim([min(X),max(X)])
    axis.set_xlabel(Xlabel,fontsize=FS)
    axis.set_ylabel(Ylabel,fontsize=FS)
    #axis.legend(handles = [el.params[3]],numpoints=1,loc=2)
    axis.set_xticks(xticks)
    axis.vlines(xticks, ymin = min(Y), ymax = max(Y))
    
    if plotA == 'plot':
        axis.plot(X,Y,el.params[1],linewidth=lw)
    if plotA == 'scatter':
        axis.plot(X,Y,el.params[2],markersize = ms)
 
    return axis

############################################################

# LAB FUNCTIONS

############################################################

def Eq1(Z,m,n): # Boar energy
    # Z atomic number
    # from m to n orbital
    # negative sign means photon emitted
    # positive sign means photon absorbed
    return R * Z**2 * ( (1/m**2) - (1/n**2) ) * (1/1000) # keV

def Eq2(Z,C,sigma): # Mosley's law
    # Z atomic number
    # C and sigma are constants for a given series (m,n pair)
    return C * (Z-sigma)**2 # eV?

def Eq3(Z,n): # spin orbit energy
    # Z atomic number
    return (Z**n)/1E7 # keV?

def Eq5(E,n): # power law describing mass absorption coefficient
    # n is an as yet un-known power
    return E**(-n) # (cm^2 g^-1 ?)

def Eq4(case,x,lam,E,n): # the attenuation of X-rays as they pass through matter
    # x is the distnace traveled through material (cm)
    # rho is density of material passed through (g cm^-3)
    #rho_Zr = 6.44 #g cm^3
    # mu(E_gamma) is the absoption coefficient, strongly energy dependant (cm^2 g^-1 ?)
    # mu goes as E^{-n} - Eq 5
    I0 = Mo_I0 # CPS, found using the CPS corresponding to Mo's K_alpha in PartIII
    if case == 1:
        return I0 * np.exp(-x/lam) # CPS?
    if case == 2:
        return I0 * np.exp(-Eq5(E,n) * rho_Zr * x)

def Edge(path):
    E = []

    LL,UL = 4,2052
    
    for el in Samples:   
        data = OpenCSV(path+el.title,LL,UL,int)
        X,Y = XY(True,data,0,1)
    
        alpha = [ el.alpha , Y[(np.abs(el.alpha - X)).argmin()] ]
        beta = [ el.beta , Y[(np.abs(el.beta - X)).argmin() ] ]
        E.append(alpha)
        E.append(beta)
        
        if el.doublet != 'no':
            doub = [ el.doublet[0] , Y[(np.abs(el.doublet[0] - X)).argmin()] ]
            E.append(doub)
    
    X2,Y2 = [],[]
    for e in E:
        X2.append(e[0])
        Y2.append(e[1])
    
    return X2,Y2
    
        


    
########################################################################
    
# FITTING FUNCTIONS

########################################################################

def Mosley(el,N):
    try_C = np.linspace(C_alpha*.8,C_alpha*1.2,N)
    try_sigma = np.linspace(sigma_alpha*.5,sigma_alpha*1.5,N)
    
    pos = []
    for c in try_C:
        for s in try_sigma:
            pos.append([c,s])
    
    def minfunc(pos):
        c,s = pos
        return np.abs(el.alpha - Eq2(el.Z,c,s))
    
    c,s = min(pos, key=minfunc)
    z = minfunc([c,s])
       
    return c,s,z
mos_C = np.array([0.0102272329008, 0.0103579845799, 0.0103130386903, 0.0104111024496, 0.0106971217477, 0.0109095932262]).astype(float)
mos_Sigma = np.array([0.944444444444, 1.04554554555, 0.83033033033, 0.864364364364, 1.14364364364, 1.15465465465]).astype(float)
mos_error = np.array([2.17700431904e-07, 2.77797467518e-06, 6.6285794631e-07, 9.65495324579e-07, 9.08419715984e-06, 9.28880223938e-06]).astype(float)

def DoubletFit(el):    
    E_so = np.abs(el.beta - el.doublet[0])    
    ns = np.linspace(3,5,1000)
    
    def minfunc(n):
        return np.abs(E_so - Eq3(el.Z,n))
    
    n = min(ns, key=minfunc)
    err_E = minfunc(n)
    err_E_per = (err_E/E_so) * 100        
    result = [n,err_E_per]    
    return result
    
def LamFit(I,x):
    lams = np.linspace(0,10,1000)
    I0 = Mo_I0
    
    def minfunc(l):
        return np.abs((I0-I) - Eq4(1,x,l,0,0))
    l = min(lams, key=minfunc)
    
    l1 = min(lams, key=minfunc)
    lams2 = np.linspace(l1*.8,l1*1.2,1000)
    
    l2 = min(lams2, key=minfunc)
    lams3 = np.linspace(l2*.9,l2*1.1,1000)
    
    l = min(lams3, key=minfunc)
    Att = Eq4(1,x,l,0,0)
    y = Mo_I0 - Att
    err_y = np.abs(I-y)
    err_y_per = (err_y/y)*100
    
    return [l,Att,y,err_y,err_y_per]

    
        
###################################################################
    
# ERROR FUNCTIONS
    
#####################################################################

def d_keV(keV): # error on the x-axis in keV
    channel = KeVToBin(A,B,keV)
    return dA + dB*channel

def d_mu(t):
    lam,err_lam = t.lam[0],t.lam[3]
    return (err_lam/rho_Zr) * (1/lam**2)
    
    
