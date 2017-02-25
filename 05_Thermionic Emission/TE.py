from numpy import polyfit,poly1d,linspace,std,array,float64,pi,sqrt,logspace,dot,log10,exp,average
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import csv
from scipy import optimize

# Tables
TT = [1000,1100,1200,1300,1400,1500,
      1600,1700,1800,1900,2000,2100,2200,2300,2400,
      2500,2600,2700,2800,2900,3000,3100,3200,3300,
      3400,3500,3600] # True Temp (K)
      
TE = [.114,.128,.143,.158,.175,
      .192,.207,.222,.236,.249,.260,.270,.279,.288,
      .296,.303,.311,.318,.323,.329,.334,.337,.341,
      .344,.348,.351,.354] # Total Emissivity

      
BT = [966,1059,1151,1242,1332,
      1422,1511,1599,1687,1774,1861,1946,2031,2115,
      2198,2280,2362,2443,2523,2602,2681,2759,2837,
      2913,2989,3063,3137] # Brightness Temp (K)
      
TABLE = []
for i in range(27):
    row = [TT[i],TE[i],BT[i]]
    TABLE.append(row)
ep_av = average(TE[9:15])

XBT_1 = [] # Experimental Brightness Values in K
XBT = [] # Experimental Brightness Values in K
with open('ExpBTemp.csv', newline='',encoding='utf-8') as d:
    reader = csv.reader(d)
    for x in reader:
        XBT_1.append(x)
for x in XBT_1:
    XBT.append(float(x[0]))

# Functions
def Ex2Data(num): # imports data from the second experiment
    DATA = []
    name = 'Thermionic Emission - %s.csv' % num
    with open(name, newline='', encoding='utf-8') as d:
        reader = csv.reader(d)
        for row in reader:
            DATA.append(row)
    for row in DATA[2:]:
        row[0] = float(row[0]) # V_anode (V)
        row[1] = float(row[1]) # I_fil (A)
        row[2] = float(row[2]) # I_fil (mA)
    for row in DATA[:1]:
        row[0] = float(row[0])
    return DATA # Values in (V,mA)
I159 = Ex2Data(1.59)
I178 = Ex2Data(1.78)
I200 = Ex2Data('2.00')
I212 = Ex2Data(2.12)
I220 = Ex2Data(2.2)
I229 = Ex2Data(2.29)

Is = [I159,I178,I200,I212,I220,I229]
    
def Leg(C,M,L): # use for legend on graph
    return mlines.Line2D([],[],color=C,marker=M,markersize=15,label=L)

def best_fit_Find(X, Y):
    
    xbar = sum(X)/len(X)
    ybar = sum(Y)/len(Y)
    n = len(X) # or len(Y)

    numer = sum([xi*yi for xi,yi in zip(X, Y)]) - n * xbar * ybar
    denum = sum([xi**2 for xi in X]) - n * xbar**2

    b = numer / denum
    a = ybar - b * xbar

    print('best fit line:\ny = {:.2f} + {:.2f}x'.format(a, b))

    return a, b
# results from best_fit_Find(BB,TT)
def BestFit(x):
    return (1.2 * x) - 201.68

def BestFitGraph():
    LL = 950
    UL = 3150

    X = linspace(LL,UL,100)
    
    plt.figure(1,figsize = (30,15))
    plt.xlabel("Brightness Temp (K)",fontsize = 20)
    plt.ylabel("True Temp (K)",fontsize = 20)
    plt.gca().set_xlim(left = LL, right = UL)
    plt.grid(True)
    
    TD = Leg('b','s','Table Data')
    BF_p = Leg('r',"",'Python Best Fit')
    BF_l = Leg('g',"",'$(1.2 * T_B(K)) - 201.68(K)$')
    plt.legend(handles=[TD,BF_p,BF_l], loc = 2, numpoints = 1)
    
    Z = polyfit(BT,TT,3)
    y2 = poly1d(Z)
    plt.plot(BT,TT,'bs', X,y2(X),'r--', X,BestFit(X),'g--')
    return
 
def Error(arg):
    LL = 950
    UL = 3150
    
    ERROR = []
    X = linspace(LL,UL,100)
    Z = polyfit(BT,TT,3)
    
    y2 = poly1d(Z)
    for x in X:
        row = y2(x)-BestFit(x)
        ERROR.append(row)
    ERROR = array(ERROR)
    y_err = std(ERROR,dtype=float64) # TTemp_err
    x_err = std(array([1902.1,1848.1,1920.1])) # BTemp_err
    
    if arg == True:
        print("Error in Brightness Temp = %s K" % x_err)
        print("Error in True Temp = %s K" % y_err)
        
    return x_err,y_err

# Error Variables and Functions
BTemp_err,TTemp_err = Error(False)
Amp_err = .005
Volt_err = .005
V0_err = .1
I0_err = [1E-5,1E-4,1E-2,1E-3,1E-3]

def Power_err(I,V): # I and V are arrays of same length
    ERR = []
    for i in range(len(I)):
        err = (float(I[i]) * Volt_err) + (Amp_err * float(V[i])) + (Amp_err * Volt_err)
        ERR.append(err)
    ERR = array(ERR)
    return ERR

def J0_T2_err(dI,dT,I,T): # I and T are arrays with length of T = length of I+1
    ERR = []
    dt = dT
    Tmag = sqrt(dot(T,T))
    Imag = sqrt(dot(I,I))
    for i in range(len(dI)):
        da = float(dI[i])
        a = float(I[i])
        t = float(T[i])
        err = (a+da)/(t**2 + 2*t*dt + dt**2)
        ERR.append(err)
    ERR = array(ERR)
    return ERR * ((Tmag**2)/Imag)

def W_err(J0,T): # and J0 and T are arrays
    J0 = array(J0)
    T = array(T)
    J0mag = sqrt(dot(I0_err,I0_err))
    Tmag = sqrt(dot(I0_err,I0_err))
    option1 = J0mag/Tmag
    option2 = std(J0/T)
    return option1
        
def TrueTempGraph(arg):
    X = array(XBT)
    Y = BestFit(X)
    
    plt.figure(2,figsize = (30,15))
    plt.errorbar(X,Y, xerr=BTemp_err, yerr=TTemp_err, fmt='ko')
    plt.xlabel("Brightness Temp (K)",fontsize = 20)
    plt.ylabel("True Temp (K)",fontsize = 20)
    plt.grid(True)
    
    if arg == True:
        print("")
        print("Error in Fig 2:")
        print("X err = %s K" % BTemp_err)
        print("Y err = %s K" % TTemp_err)
        print("")
    return

def Part1A(arg):
    if arg == True:
        print(BestFit(array(XBT)))
    return BestFit(array(XBT))
    
XTT = Part1A(False)

def Part1B(arg1,arg2):
    DATA_1 = [] # Power (W) , B Temp (K) , T Temp (K)
    X1 = [] # B Tem
    X2 = [] # T Tem
    Y = [] # Power
    I = []
    V = []
    
    with open('Thermionic Emission - PowerTemp.csv', newline='',encoding='utf-8') as d:
        reader = csv.reader(d)
        for row in reader:
            DATA_1.append(row)
    for row in DATA_1[1:]:
        X1.append(float(row[1]))
        X2.append(float(row[2]))
        Y.append(float(row[0]))
        V.append(float(row[4]))
        I.append(float(row[3]))
    
    A = 1.30E-05 # area m^2 
    sigma = 5.68E-08 # W m^-2 K^-4
    
    def SB_1(eps,T):
        power = eps * A * sigma * T**4
        return power
    
    def SB_2(ep,n,T):
        power = ep * A * sigma * T**n
        return power
    
    def FindFit1(X,Y):  
        alpha = A * sigma
        def fitfunc(p, x):
            return p[0] * alpha * (x**4)
        def errfunc(p, x, y):
            return y - fitfunc(p, x)

        xdata=array(X)
        ydata=array(Y)

        qout,success = optimize.leastsq(errfunc, [.3], args=(xdata, ydata),maxfev=3000)

        out = qout[:]
        out[0] = qout[0]
        return out[0]
    ep = FindFit1(X1,Y)
    
    def Num1():
        LL = 1000
        UL = 2200
        X3 = linspace(LL,UL,100)
        Y2 = SB_1(ep,X3)
            
        plt.figure(3,figsize = (30,15))
        Dleg = Leg('b','s','Data')
        Y2Leg = Leg('r',"",'$\epsilon = %s$' % round(ep,2))
    
        plt.legend(handles=[Dleg,Y2Leg], loc=2, numpoints = 1)
        plt.errorbar(X1,Y, xerr=BTemp_err, yerr=Power_err(I,V), fmt='bs', label='Data')
        plt.plot(X3,Y2,'r-')        
        plt.xlabel("Brightness Temp (K)",fontsize = 20)
        plt.ylabel("Power (W)",fontsize = 20)
        plt.grid(True)
        
        if arg1 == True:
            print("Results for Fig 3:")
            print("Data Fit ep = %s" % ep)
            print("Known Average ep for range = %s" % ep_av)
            print("epsilon error = %s" % (ep_av-ep))
            print("")
            print("Brightness Temp Error = %s K" % BTemp_err)
            print("Power Error = %s W" % Power_err(I,V))
            print("")
        return
    
    def Num2():     
        def FindFit2(X,Y):
            alpha = ep * A * sigma
            def fitfunc(p, x):
                return alpha * (x ** p[0])
            def errfunc(p, x, y):
                return y - fitfunc(p, x)

            xdata=array(X)
            ydata=array(Y)

            qout,success = optimize.leastsq(errfunc, [4],
                               args=(xdata, ydata),maxfev=5000)

            out = qout[:]
            out[0] = qout[0]
            return out[0]  
        n = FindFit2(X2,Y)
            
        LL = 1000
        UL = 2200
        X3 = linspace(LL,UL,100)
        Y2 = SB_2(ep,n,X3)
        Y3 = (ep*A*sigma)* X3**4
        
        plt.figure(4,figsize = (30,15))
        Dleg = Leg('b','s','Data')
        Y2Leg = Leg('r',"",'Fit: $\epsilon A \sigma T^{%s}$' % round(n,3))
        Y3Leg = Leg('g',"",'$\epsilon A \sigma T^4$')
        plt.legend(handles=[Dleg,Y2Leg,Y3Leg], loc=2, numpoints = 1)
        plt.errorbar(X1,Y, xerr=TTemp_err, yerr=Power_err(I,V), fmt='bs')
        plt.plot(X3,Y2,'r-', X3,Y3,'g-')        
        plt.xlabel("True Temp (K)")
        plt.ylabel("Power (W)")
        plt.grid(True)
        
        if arg2 == True:
            print("Results for Fig 4:")
            print("Experimental n = %s" % n)
            print("Error in n = %s" % (4-n))
            print("True Temp err = %s K" % TTemp_err)
            print("Power err = %s W" % Power_err(I,V))
            print("")
            
        return
    return Num1(),Num2()

def Part2(arg1,arg2,arg3):
    l = .0317 # length of filiment - m
    beta_n2 = .93 # num
    #eta = 1.758820024E11 # e/m - C/kg
    diam = 1.58 # diameter of anode - cm
    r_a = (diam/2)/100 # radius of anode - m
    ep_0 = 8.854E-12 # F/m
    eta_goal = (1.602E-19)/(9.101E-31) # MKS
    
    alpha = (8 * pi * ep_0 * l * sqrt(2) * beta_n2)/(9 * r_a)
    conversion = (1/sqrt(5.61E35)) * (2.998E8) * (1/sqrt(1.602E-19))
    
    def Child(x,eta,n): #
        return alpha * conversion * sqrt(eta) * x**n
    
    def FindFit(X,Y): # e/m * V_0^n
            
        def fitfunc(p, x): #
            return (alpha * conversion) * sqrt(p[0]) * x**p[1]
        def errfunc(p, x, y):
            return y - fitfunc(p, x)

        xdata=array(X)
        ydata=array(Y)

        qout,success = optimize.leastsq(errfunc, [1.7E11,1.5],args=(xdata, ydata),maxfev=5000)

        out = qout[:]
        out[0] = qout[0]
        out[1] = qout[1]
        return [out[0],out[1]]
            
    def XY(Data):
        X = []
        Y = []
        for row in Data[3:]:
            X.append(float(row[0]))
            Y.append(float(row[1])) # [2] for mA, [1] for A
        return X,Y
    X0,Y0 = XY(I159)
    X1,Y1 = XY(I178)
    X2,Y2 = XY(I200)
    X3,Y3 = XY(I212)
    X4,Y4 = XY(I220)
    X5,Y5 = XY(I229)
    
    ymin = 1E-4
    xmin = .1
    
    def plot(num):
    
        L0 = Leg('r','o','$I_{plate} = 1.59\ mA$')
        L1 = Leg('y','o','$I_{plate} = 1.78\ mA$')
        L2 = Leg('g','o','$I_{plate} = 2.00\ mA$')
        L3 = Leg('c','o','$I_{plate} = 2.12\ mA$')
        L4 = Leg('b','o','$I_{plate} = 2.20\ mA$')
        L5 = Leg('m','o','$I_{plate} = 2.29\ mA$')
    
        plt.figure(num,figsize = (15,15))
        plt.legend(handles=[L0,L1,L2,L3,L4,L5], loc=2, numpoints = 1)
        plt.plot(X0,Y0,'ro', X1,Y1,'yo', X2,Y2,'go', X3,Y3,'co', X4,Y4,'bo', X5,Y5,'mo')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel("$V_{anode}\ (V)$",fontsize = 24)
        plt.ylabel("$I_{plate}\ (A)$",fontsize = 24)
        plt.gca().set_ylim(bottom = ymin)
        plt.gca().set_xlim(left = xmin)
        plt.grid(True)
        return
     
    def Num1():
        
        if arg1 == True: plot(5)
        
        eta2,n2 = FindFit(X2[:-4],Y2[:-4])
        eta3,n3 = FindFit(X3[:-3],Y3[:-3])
        eta4,n4 = FindFit(X4[:-2],Y4[:-2])
        eta5,n5 = FindFit(X5[:-1],Y5[:-1])
        
        X = logspace(-2,3,100)
        Y_2 = Child(X,eta2,n2)
        Y_3 = Child(X,eta3,n3)
        Y_4 = Child(X,eta4,n4)
        Y_5 = Child(X,eta5,n5)
        V_0 = 1.2
        if arg1 == True: 
            plt.plot(X,Y_5,'m-', X,Y_4,'b-', X,Y_3,'c-', X,Y_2,'g-', V_0,ymin,'mo',markersize = 20)
            print("Results from Fig 5:")
            print("Extraplated value for V_0 = %s V" % V_0)
            print("V_0 err = %s V" % V0_err)
            print("Voltage err = %s V" % Volt_err)
            print("Current err = %s A" % Amp_err)
            print("error not shown on graph because it looks sloppy")
            print("")
    
        return eta5,n5,V_0
        
    eta,n,V_0 = Num1()
    eta_err = eta_goal-eta
    n_err = (3/2)-n
    if arg1 == True:
        print("e/m and n error:")
        print("em = %s E11 C/kg" % (eta/1E11))
        print("Error in e/m = %s E11 C/kg" % (eta_err/1E11))
        print("Error in e/m = %s percent" % ((eta_err/eta_goal)*100))
        print("n = %s" % n)
        print("Error in n = %s" % n_err)
        print("")
    
    def Num2():
        
        if arg2 == True: plot(6)
        
        eta0,n0 = FindFit(X0,Y0)
        eta1,n1 = FindFit(X1,Y1)
        eta2,n2 = FindFit(X2[1:],Y2[1:])
        eta3,n3 = FindFit(X3[3:],Y3[3:])
        eta4,n4 = FindFit(X4[4:],Y4[4:])
        # eta5,n5 does not show a horizontal asymptote
        
        X = logspace(-2,3,100)
        Y_0 = Child(X,eta0,n0)
        Y_1 = Child(X,eta1,n1)
        Y_2 = Child(X,eta2,n2)
        Y_3 = Child(X,eta3,n3)
        Y_4 = Child(X,eta4,n4)
        
        I0 = [.000135,.00082,.004,.0089,.0122]
        
        if arg2 == True: 
            plt.plot(X,Y_0,'r-', X,Y_1,'y-', X,Y_2,'g-', X,Y_3,'c-', X,Y_4,'b-')
            print("Results from Fig 6:")
            print("I_0 = %s mA" % (array(I0)*1E3))
            print("I_0 uncertainty = %s mA" % (array(I0_err)*(1E3)))
            print("Voltage err = %s V" % Volt_err)
            print("Current err = %s A" % Amp_err)
            print("error not shown on graph because it looks sloppy")
            print("")
        
        MS = 20
        
        if arg2 == True: 
            plt.plot(xmin,I0[0],'ro', xmin,I0[1],'yo', xmin,I0[2],'go', xmin,I0[3],'co', xmin,I0[4],'bo', markersize=MS)
        
        I0 = array(I0)
        A = 2 * pi * r_a * l
        J0 = I0/A
        
        return J0
    
    def Num3():
        
        J0 = Num2()
        T_data = []
        for I in Is[:-1]:
            for row in I[1]:
                if len(row) > 1:
                    T_data.append(float(row)) 
        T_data = array(T_data)
        
        def XY():
            X = []
            Y = []
            for t in T_data:
                x = 1000/float(t)
                X.append(x)
            for i in range(5):
                y = (J0[i]/(X[i])**2)*10**(-3)
                Y.append(y)  
            return X,Y
        X_data,Y_data =  XY()
        
        def Work_Function(X,Y,n): # X and Y are equal length arrays, n is the final index
            X = array(X)
            Y = array(Y)
            
            k =  1.3807E-23 # J/K
            e = 1.60217662E-19 # C
            #A_0 = 1.2E6 # A m^-2 K^-2
            
            yf,yi = log10(Y[n]), log10(Y[0])
            xf,xi = X[n],X[0]
            
            m = abs(1000 * ((yf-yi)/(xf-xi))) # slope # K
            W = m * (k/e) * (e/1.602E-19) # V, e's cancel out when converting Joule to eV
            return m,W
        
        m,W = Work_Function(X_data,Y_data,4)
        
        x_err = 1/TTemp_err**2
        y_err = J0_T2_err(I0_err,TTemp_err,J0,T_data)
        m_err = std(log10(y_err)/x_err)
        
        def plot(num):
            
            Dleg = Leg('b','s','Data')
            Fleg = Leg('r',"",'Linear Fit')
            
            plt.figure(num,figsize = (8,15))
            plt.legend(handles=[Dleg], loc=1, numpoints = 1)
            # y_err is not put on graph because it looks cluncky in semi-log graph
            plt.errorbar(X_data,Y_data, xerr=x_err, fmt='bs')
            plt.yscale('log')
            plt.xlabel("1000/T",fontsize = 20)
            plt.ylabel("$J_0/T^2$ per unit area $[A\ K^{-2}\ m^{-2}]$ x $10^{-4}$",fontsize = 20)
            plt.grid(True)
            return
        
        if arg3 == True: 
            plot(7)
            print("Error in Fig 7:")
            print("Y err not shown on graph; looks cluncky in semi-log plot")
            print("Error in Temp = %s K" % TTemp_err)
            print("Error in J_0/T^2 = %s A K^-2 m^-2$" % y_err)
            print("")
            print("Slope = %s K" % m)
            print("Slope error = %s K" % m_err)
            print("Work Function for Tungsten = %s V" % W)
            print("Error in Work Function = %s V" % W_err(J0,T_data))
            print("")
    
    return Num1(),Num2(),Num3()  

def Call():
    BestFitGraph()
    Error(True)
    TrueTempGraph(True)
    Part1A(True)
    Part1B(True,True)
    Part2(True,True,True)
    return
Call()