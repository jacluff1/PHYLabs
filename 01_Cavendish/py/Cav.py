import numpy as np
import matplotlib.pyplot as plt 
plt.style.use('~/.config/matplotlib/default.mplstyle')
import djak.gen as dg
import djak.plot as dp
import djak.fits as df

# Classes----------------------------------------------------------

class ball:
    def __init__(self,m,d,conv): # g -> kg , cm -> m
        self.conv = conv
        if self.conv == True:
            self.m = m/1000 # g -> kg
            self.d = d/100 # cm -> m
        else:
            self.m = m # kg
            self.d = d # m
        self.r = self.d/2 # m

# error functions

def G_s_err(k,theta_s,R,M,m,d):
    err1 = (k.err/k.val)**2
    err2 = (theta_s.err/theta_s.av)**2
    err3 = 2 * (R.err/R.val)**2
    err4 = (M.val/M.err)**2
    err5 = (m.val/m.err)**2
    err6 = (d.val/d.err)**2
    return np.sqrt(err1 + err2 + err3 + err4 + err5 + err6)

def G_d_err(k,theta_d,R,M,m,d):
    err1 = (k.err/k.val)**2
    err2 = (theta_d.err/theta_d.val)**2
    err3 = 2 * (R.err/R.val)**2
    err4 = (M.val/M.err)**2
    err5 = (m.val/m.err)**2
    err6 = (d.val/d.err)**2
    return np.sqrt(err1 + err2 + err3 + err4 + err5 + err6)

def k_err(T,I):
    err1 = 2 * (T.av/T.err)**2
    err2 = (I.err/I.val)**2
    return np.sqrt(err1 + err2)

def I_err(m_m, d, m_r, m_b, l_b, w_b):
    eta = d.val**2 + (2/5)*m_r.val**2
    d_eta = np.sqrt( 2*d.err**2 + 2*m_r.err**2 )
    d_Is = np.sqrt( (m_m.err/m_m.val)**2 + (d_eta/eta)**2 )

    mu = l_b.val**2 + w_b.val**2
    d_mu = np.sqrt( 2*l_b.err**2 + 2*w_b.err**2 )
    d_Ib = np.sqrt( (m_b.err/m_b.val)**2 + (d_mu/mu)**2 )

    d_I = np.sqrt( d_Is**2 + d_Ib**2 )
    return d_I

def R_err(W,w_t,M_d):
    return np.sqrt(W.err*2 + w_t.err**2 + M_d.err**2)

def d_err(l_b,m_d):
    return np.sqrt(l_b.err**2 + m_d.err**2)

def beta_err(gamma,T):
    err1 = (gamma.err/gamma.val)**2
    err2 = (T.err/T.av)**2
    return np.sqrt(err1 + err2)

def theta_laser_err(theta,theta_err,C):
    return np.sqrt( (theta_err/theta)**2 + (C.err/C.val)**2 )
        
# GLOBAL VARIABLES-------------------------------------------------

nppath = "../npy/"
csvp = "../csv/"
gp = "../graphs/"
per = '$\%$'

csvs = ["angle_cal","static1","static2","driven","equilib","extra_decay"]
npys = csvs

M1 = ball(1054.5,5.15,True)
M2 = ball(1052.9,5.15,True)
m1 = ball(14.5,.85,True)
m2 = ball(14.5,.85,True)
M = ball((M1.m+M2.m)/2,(M1.d+M2.d)/2,False)
m = ball((m1.m+m2.m)/2,(m1.d+m2.d)/2,False)

S = dg.var('S',np.array([43.6,28.4,18.3,11.9,7.8,5.1])/100, .5/100, 'meters')
L = dg.var('L',762/100, .1, 'meters')
T = dg.var('T',[396., 385., 380., 380., 382., 382.5], 'std', 's') # Period

theta_s = dg.var('theta_s',np.array([1.88152248, 1.70676281, 1.73349876, 1.78764081])/1000, 'std', 'rad') 
theta_d = dg.var('theta_d',.001,.007,'mRad')

G_s = dg.var('G_s', np.array([5.64897607e-11, 5.12428757e-11, 5.20455807e-11, 5.36711108e-11]), 52684376333.6, 'N m^2 kg^-2' )
G_d = dg.var('G_d', 3.01479199888e-11, 52684376333.6, 'N m^2 kg^-2')


# boom
l_b = dg.var('l_b',14.9/100, .1/100, 'm') # length of boom
m_b = dg.var('m_b',7.1/1000, .5/1000, 'kg') # mass of boom
w_b = dg.var('w_b',.1/100, .05/100, 'm') # width of boom

# cavendish box
w_t = dg.var('w_t',.1/100, .05/100, 'm') # glass thickness
W = dg.var('W',2.5/100, .05/100, 'm') # width of box

# m and m'
m_m = dg.var('m_m',m.m, .5/1000, 'kg') # mass of m
m_d = dg.var('m_d',m.d, .05/100, 'm') # diameter of m
m_r = dg.var('m_r',m_d.val/2, m_d.err/2, 'm') # radius of m

# M and M'
M_m = dg.var('M_m',M.m, .5/1000, 'kg') # mass of M
M_d = dg.var('M_d',M.d, .05/100, 'm') # diameter of M
M_r = dg.var('M_r',M_d.val/2, M_d.err/2, 'm') # radius of M

# fitting params
#theta_0, A, gamma, omega
par_m = [0,10,1.5,(2*np.pi)/(.385)]
theta_e = dg.var('theta_e',0.0, 0.1, 'mRad') # equilibrium angle
A = dg.var('A',10.5, 0.2, 'mRad') # initial amplitude of decaying wave
gamma = dg.var('gamma',1.6/1000, 0.1/1000, 's^-1')
omega = dg.var('omega',16.3/1000, 0.2/1000, 's^-1') # angular frequency 
# misc
d = dg.var('d',(l_b.val-m_d.val)/2, d_err(l_b,m_d), 'm') # distance from rotation axis to center of small sphere
I_b = (m_b.val/12) * (l_b.val**2 + w_b.val**2)
I_s = (2*m_m.val)*(d.val**2 + (2/5)*m_r.val**2)
I = dg.var('I',I_b + I_s, I_err(m_m,d,m_r,m_b,l_b,w_b), 'kg m^2') # moment of inertia of boom, m and m'
beta = dg.var('beta',2*T.av*gamma.val, beta_err(gamma,T), None) # decay parameter - 
R = dg.var('R',W.val/2 + w_t.val + M_d.val/2, R_err(W,w_t,M_d), 'm') # distance between center of Ms and m's
C = dg.var('C',.89,.01,None) # calibration constant

# DATA ------------------------------------------------------------

def Write(): # writes all .csv files to .npy files
    for f in csvs:
        dg.CSV(csvp,nppath,f)

def AuxOne(room_L,laser_S, theta_comp_ext): # returns coordinates for laser/comp calibration plot
    # room_L is the distance between boom mirror and opposite wall
    # laser_S is the distance between adjacent extrema shown by laser on opposite wall
    # theta_comp_extrema is the theta_comp of each extrema
    theta_comp = []
    for i in range(len(theta_comp_ext[1:])):
        theta_comp.append(theta_comp_ext[i+1] - theta_comp_ext[i])
    theta_comp = abs(np.array(theta_comp))

    def Laser(s):
        result1 = (1/2) * np.arctan(s/L.val)
        return result1*1000 
    theta_laser = Laser(S.val)

    return theta_comp,theta_laser 

def AuxTwo(t1,t2,type_t,X_data,Y_data,dx): # approximate times of applicable static measurement (turning points on either end)
    # t1 is approximate time of first turning point
    # t2 is approximate time of the last turning point
    # type_t specifies t1 and t2 as either maxima or dg.Minima
    # X_data is the time array of data
    # Y_data is the theta_laser array of data
    # dx is the range to test for true turning point
    if type_t == 'maxima':
        imin,xmin,ymin = dg.maxima(X_data,Y_data,t1,dx,ret='all')
        imax,xmax,ymax = dg.maxima(X_data,Y_data,t2,dx,ret='all')
    if type_t == 'minima':
        imin,xmin,ymin = dg.minima(X_data,Y_data,t1,dx,ret='all')
        imax,xmax,ymax = dg.minima(X_data,Y_data,t2,dx,ret='all')
    av_y = np.average(Y_data[imin:imax])
    Y_av = []
    for i in range(len(X_data[imin:imax])):
        Y_av.append(av_y)
    return xmin,ymin,xmax,ymax,np.array(X_data[imin:imax]),np.array(Y_av)
    
    
# LAB EQUATIONS----------------------------------------------------

def Theta_s(theta_1,theta_2): # returns static displacement angle
    # theta_1 = displacement angle 1
    # theta_2 = displacement angle 2
    theta_1 = theta_1[0]
    theta_2 = theta_2[0]
    return (1/2)*abs(theta_2-theta_1)

def Period(t1,t2,conv,n): # returns period
    # t1 = beginning of period 1
    # t2 = end of period n
    # conv = conversion factor (ex ks -> s, conv = 1000)
    # n = number of periods
    return abs(t2-t1)*(conv/n)
    
def K_s(I,T):
    # I = moment of inertia
    # T = period
    return (2*np.pi/T)**2 * I
k = dg.var('k',K_s(I.val,T.av), k_err(T,I), 'N m') # torsion constant

def G(theta): # return gravitational constant
    #print(type(k.val),type(theta),type(R.val),type(M_m.val),type(m_m.val),type(d.val))
    return ( k.val * theta * R.val**2 )/( 2 * M_m.val * m_m.val * d.val)
def G2(theta):
    return (d.val**2 *omega.val**2 * l_b.val * theta)/(2*M_m.val)

def dampedWave(X,Y,par_m):
    # par = [theta_0, A, gamma, omega]
 
    def envelope(x):
        t = x-X[0]
        return par_m[1] * np.exp(-par_m[2]*t)
    
    def fitfunc(p,x):
        t = x-X[0]
        return p[0] + p[1]*np.exp(-p[2]*t)*np.cos(p[3]*t)
    
    par_f = df.least_square(fitfunc,par_m,X,Y)
    #par_f[3] += .2
        
    Y_env = envelope(X)
    Y_man = fitfunc(par_m,X)
    Y_fit = fitfunc(par_f,X)
    
    return Y_env,Y_man,Y_fit,par_f

def theta_driven(X,Y,theta_1,printA=False):
    x1 = np.e**(-(gamma.val*T.av)/2)
    
    t_max_approx = np.array([3.26, 3.64, 4.04, 4.41, 4.79, 5.17]) # 2.87
    t_min_approx = np.array([3.06, 3.44, 3.82, 4.20, 4.59, 4.98]) # 5.36
    i_max,i_min = [0],[]
    for t in t_max_approx:
        i_max.append(dg.maxima(X,Y,t,25))
    for t in t_min_approx:
        i_min.append(dg.minima(X,Y,t,25))
    i_ext = np.hstack((np.array(i_max),np.array(i_min)))
    i_ext = np.sort(i_ext)
    THETA = Y[i_ext]
    N = len(THETA) 
    
    # needs work
    def theta_n(n,x):
        return theta_e.val + (theta_1 - theta_e.val)*(-x)**(n-1)
    THETA2 = []
    for n in np.arange(1,N+1):
        THETA2.append(theta_n(n,x1))
    THETA2 = np.array(THETA2)
    
    Xs = []
    for i in np.arange(1,N-2):
        Xs.append( -(THETA[i+2] - THETA[i+1])/(THETA[i+1] - THETA[i]) )
    x2 = dg.var('x2',Xs,'std','mRad')
    
    x3_val = 1 - (THETA[1] - THETA[11])/(THETA[1]-THETA[2]+THETA[3]-THETA[4]+THETA[5]-THETA[6]+THETA[7]-THETA[8]+THETA[9]-THETA[10])
    def x3_err(theta_err):
        return theta_err * (1-x3_val) * np.sqrt( (13-1)*(1-x3_val)**2 + 2*x3_val ) / abs(THETA[0]-THETA[12])
    x3 = dg.var('x3',x3_val,np.average(x3_err(theta_laser_err(THETA,.05,C))),None)
    
    theta_d_val = ( (1-x3.val)*(THETA[0]-THETA[1]+THETA[2]-THETA[3]+THETA[4]-THETA[5]+THETA[6]-THETA[7]+THETA[8]-THETA[9]+THETA[10]-THETA[11]+THETA[12]) - THETA[0] + THETA[12] ) / ( (13-1) * (1+x3.val) )
    
    err1 = np.average(theta_laser_err(THETA,.05,C)) * np.sqrt( (13-1)*(1-x3.val)**2 + 2*x3.val ) / ( (13-1)*(1+x3.val) )
    err2 = x3.err * ( 2 *(THETA[0]-THETA[1]+THETA[2]-THETA[3]+THETA[4]-THETA[5]+THETA[6]-THETA[7]+THETA[8]-THETA[9]+THETA[10]-THETA[11]) + (THETA[12] - THETA[1] ) ) / ( (13-1)*(1+x3.val)**2 )
    theta_d_err = np.sqrt( err1**2 + err2**2 )
    theta_d = dg.var('theta_d',theta_d_val,theta_d_err,'mRad')
    
    if printA == True:
        print(theta_d.name,theta_d.val,theta_d.err,theta_d.units)
    
    return theta_d
    
# Graph Functions--------------------------------------------------

def Angle_calibration(plot1=False,plot2=False,printA=False):
    plt.close('all')
    angle_data = np.load(nppath+"angle_cal.npy")
    X = angle_data[:,0] # sec
    Y = angle_data[:,1] # mRad

    x_minima = [62,178,293,407] # sec
    x_maxima = [120,235,350] # sec
    minima,maxima = dg.extrema(X,Y,10,x_minima,x_maxima)
    ext = minima + maxima
    ext = dg.sort(np.array(ext),0)

    if plot1 == True: # '0_1_Angles_cal', angles for angle calibration
        # plots: X , Y , style , label
        p1 = dp.data(X,Y,'-','b','$\\theta_{comp}$')
        p2 = dp.data(ext[:,0],ext[:,1],'r*','r','Extrema')
        ax = dp.ax([p1,p2],111, 'time [sec]', "$\\theta_{comp}$ [mRad]", "$\\theta_{comp}$ for Angle Calibration")
        
        dp.plot([ax],name = gp+'0_1_Angles_cal') 
        
    X_data,Y_data = AuxOne(L.val,S.val,ext[:,1]) 
    m,m_err = df.m_exp(X_data,Y_data),df.sig_m(X_data,Y_data)
    b,b_err = df.b_exp(X_data,Y_data),df.sig_b(X_data,Y_data)
    m_err_p, b_err_p = m_err*100/m, b_err*100/b
    r2 = df.R2(X_data,Y_data)   
    X_lin = np.linspace(0,max(X_data)+10,1000)
    Y_lin = df.lin_fit(X_data,Y_data,X_lin)
    
    if plot2 == True:
        p1 = dp.data(X_data,Y_data,'bs','b','Laser vs Comp')
        p2 = dp.data(X_lin,Y_lin,'r-','r','Linear Fit')
        n2 = dp.note("Y = C x + B\nm = %s $\pm$ %s\nb = %s $\pm$ %s\nR$^2$ = %s" % ( round(m,2),round(m_err,2),round(b,1),round(b_err,1),round(r2,5) ),20,7,p2.color)
        ax = dp.ax([p1,p2],111,'$\\theta_{comp}$ [mRad]','$\\theta_{laser}$ [mRad]','$\\theta_{laser}$ vs. $\\theta_{comp}$')
        ax.notes = [n2]
        
        dp.plot([ax],name = gp+'0_2_laser_vs_comp')
    
    if printA == True:
        print("")
        print("Calibrtion:")
        print("theta_comp = %s" % X_data)
        print("theta_laser = %s" % Y_data)
        print("calibration constant = %s +/- %s (%s)" % (m,m_err,m_err_p) )
        print("intercept = %s +/- %s (%s)" % (b,b_err,b_err_p) )
        print("R2 = %s" % r2 )
        print("")
    return m,m_err
#Cal,Cal_err = Angle_calibration() # multiply theta_comp by Cal
#C = dg.var(float(Cal),float(Cal_err),None) # Calibration constant
#t_err = 1 # sec

def Equilib(plot1=False,printA=False):
    D1 = np.load(nppath+'angle_cal.npy')
    D2 = np.load(nppath+'extra_decay.npy')
    X1,Y1 = D1[:,0]/100,D1[:,1]
    X2,Y2 = D2[:,0]/100,D2[:,1]

    x1,y1,x2,y2,avX,avY = AuxTwo(.55,15.76,'maxima',X2,Y2,25)
    
    d1 = dp.data(X1,Y1,'-','b','data 1')
    d2 = dp.data(X2,Y2,'-','b','data 2')
    d3 = dp.data([x1,x2],[y1,y2],'*','g','maxima')
    d4 = dp.data(avX,avY,'-','r','eq 2')
    
    if plot1 == True:
        ax1 = dp.ax([d1],121,'time [$10^2$ sec]', '$\\theta_{comp}$ [mRad]', 'Equilibrium Angle 1')
        ax2 = dp.ax([d2,d3,d4],122,'time [$10^2$ sec]', '$\\theta_{comp}$ [mRad]', 'Equilibrium Angle 2')
    
        dp.plot([ax1,ax2])
    
    eq1 = -1.8 # mRad
    eq2 = avY[0] # mRad
    
    if printA == True:
        print("eq 1 = %s mRad" % eq1)
        print("eq 2 = %s mRad" % eq2)
    
    return eq1,eq2
eq1,eq2 = Equilib()
    
def StaticG(plot=False,printA=False):
    plt.close('all')
    
    # static1
    D1 = np.load(nppath+'static1.npy')
    X1,Y1 = D1[:,0]/1000,(D1[:,1] - eq1)*C.val # time[s], theta_laser[mRad]
    
    x1_1a,y1_1a,x1_1b,y1_1b,XP_1,YP_1 = AuxTwo(.630,2.25,'maxima',X1,Y1,25)
    x1_2a,y1_2a,x1_2b,y1_2b,XP_2,YP_2 = AuxTwo(2.57,4.1,'minima',X1,Y1,25)
    x1_3a,y1_3a,x1_3b,y1_3b,XP_3,YP_3 = AuxTwo(4.40,5.92,'maxima',X1,Y1,25)
    x1_4a,y1_4a,x1_4b,y1_4b,XP_4,YP_4 = AuxTwo(6.40,7.97,'minima',X1,Y1,25)
    
    X2 = np.array([x1_1a,x1_1b,x1_2a,x1_2b,x1_3a,x1_3b,x1_4a,x1_4b])
    Y2 = np.array([y1_1a,y1_1b,y1_2a,y1_2b,y1_3a,y1_3b,y1_4a,y1_4b])
    X3 = np.hstack((XP_1,XP_2,XP_3,XP_4))
    Y3 = np.hstack((YP_1,YP_2,YP_3,YP_4))
    
    n1 = dp.note(str(round(np.average(YP_1),2)), x1_1b-.8,y1_1b+.8, 'r')
    n2 = dp.note(str(round(np.average(YP_2),2)), x1_2b-.8,y1_2b-.8, 'r')
    n3 = dp.note(str(round(np.average(YP_3),2)), x1_3b-.8,y1_3b+.8, 'r')
    n4 = dp.note(str(round(np.average(YP_4),2)), x1_4b-.8,y1_4b-.8, 'r')
    
    # static2
    D2 = np.load(nppath+'static2.npy')
    X4,Y4 = D2[:,0]/1000,(D2[:,1] - eq2)*C.val
    
    x2_1a,y2_1a,x2_1b,y2_1b,XP_5,YP_5 = AuxTwo(0.26,1.8,'minima',X4,Y4,25)
    x2_2a,y2_2a,x2_2b,y2_2b,XP_6,YP_6 = AuxTwo(2.07,3.58,'maxima',X4,Y4,25)
    
    X5 = np.array([x2_1a,x2_1b,x2_2a,x2_2b])
    Y5 = np.array([y2_1a,y2_1b,y2_2a,y2_2b])
    X6 = np.hstack((XP_5,XP_6))
    Y6 = np.hstack((YP_5,YP_6))
    
    n5 = dp.note(str(round(np.average(YP_5),2)), x2_1b-.5,y2_1b-.5, 'r')
    n6 = dp.note(str(round(np.average(YP_6),2)), x2_2b-.3,y2_2b+.3, 'r')
    
    Ts = np.array([Period(x1_1b,x1_1a,1000,4),Period(x1_2b,x1_2a,1000,4),Period(x1_3b,x1_3a,1000,4),Period(x1_4b,x1_4a,1000,4),Period(x2_1b,x2_1a,1000,4),Period(x2_2b,x2_2a,1000,4)])
    T = dg.var('T',Ts,'std','s')
    
    theta_Ds = np.array([Theta_s(YP_1,YP_2),Theta_s(YP_2,YP_3),Theta_s(YP_3,YP_4),Theta_s(YP_5,YP_6)])
    theta_s = dg.var('theta_s',theta_Ds,'std','mRad')
    
    Gs = G(theta_Ds)
    G_s = dg.var('G_s',Gs,'sts','N m^2 kg^-2')
    
    if plot == True:
        # static1
        d1 = dp.data(X1,Y1,'-','b','data 1')
        d2 = dp.data(X2,Y2,'*','g','$\\theta_D$ Period')
        d3 = dp.data(X3,Y3,'-','r','$\\theta_D$')

        ax1 = dp.ax([d1,d2,d3],121,'time [$10^3$ sec]','$\\theta_{laser}$ [mRad]','Static Measurement 1')
        ax1.notes = [n1,n2,n3,n4]

        # static2
        d4 = dp.data(X4,Y4,'-','b','data 2')
        d5 = dp.data(X5,Y5,'*','g','$\\theta_s$ Period')
        d6 = dp.data(X6,Y6,'-','r','$\\theta_s$')
        
        ax2 = dp.ax([d4,d5,d6],122,'time [$10^3$ sec]','$\\theta_{laser}$ [mRad]','Static Measurement 2')
        ax2.notes = [n5,n6]
         
        dp.plot([ax1,ax2],name=gp+'0_3_StaticG')
    if printA == True:
        print("")
        print("Static G Results:")
        print("periods = %s" % T.pval)
        print("Period: %s +/- %s %s" % (T.pav,T.perr,T.units) )
        print("Theta_s = %s %s" % theta_s.pval,theta_s.units)
        print("theta_d_static = %s +/- %s %s" % (theta_s.pav,theta_s.perr,theta_s.units))
        print("Gs = %s %s" % (G_s.val,G_s.units) )
        print("G = %s +/- %s %s" % (G_s.av,G_s.err,G_s.units) )
    return theta_s

def DrivenG(plot=False,printA=False):
    plt.close('all')
    
    D = np.load(nppath+'driven.npy')
    X,Y = D[:,0]/1000,(D[:,1] - eq2)*C.val # data
    data = dp.data(X,Y,'-','b','data') # data
    
    i_1 = dg.maxima(X,Y,2.87,25)
    x1,y1 = X[i_1],Y[i_1]
    d1 = dp.data([x1],[y1],'*','r','$X[t_0],\ Y[t_0]$')
    
    X2 = X[i_1:]
    Y2 = Y[i_1:]
    Y2_env,Y2_man,Y2_fit,par_f = dampedWave(X2,Y2,par_m)
    d2_man = dp.data(X2,Y2_man,'--','m','manual fit')
    d2_env = dp.data(X2,Y2_env,'--','c','upper envelope')
    d2_fit = dp.data(X2,Y2_fit,'-','r','least-square fit')
    
    n1 = '$\\theta_e + A e^{-\\gamma t} \cos{\\omega t}$'
    n2 = '$\\theta_e$ = %s $\pm$ %s mRad' % (round(par_f[0],1),.1)
    n3 = 'A = %s $\pm$ %s mRad' % (round(par_f[1],1),.2)
    n4 = '$\\gamma$ = %s $\pm$ %s (10$^3$ s)$^{-1}$' % (round(par_f[2],1),.1)
    n5 = '$\\omega$ = %s $\pm$ %s (10$^3$ s)$^{-1}$' % (round(par_f[3],1),.2)
    note = dp.note(n1+'\n'+n2+'\n'+n3+'\n'+n4+'\n'+n5,4,-10,'r')
    
    theta_d = theta_driven(X2,Y2,y1,printA=True)
    G_d = dg.var('G_d',G(theta_d.val),G_d_err(k,theta_d,R,M_m,m_m,d),'N m^2 kg^-2')
    dg.printvar(G_d)
    print(G_d.val,G_d.err)
    
    if plot==True:
        ax = dp.ax([data,d1,d2_env,d2_man,d2_fit],111,'time [$10^3$ sec]','$\\theta_{laser}$ [mRad]','Driven Measurement')
        ax.notes = [note]

        dp.plot([ax])#, name=gp+'0_4_DrivenG')
    
    if printA == True:
        print("Driven G Results:")
    
    return theta_d

#Angle_calibration(plot1=True)
#StaticG(plot=True)
#DrivenG()

def torun():
    #w_b = dg.var('w_b',.1/100,.05/100,'m')
    #I = dg.var('I', I_s+I_b, I_err(m_m,d,m_r,m_b,l_b,w_b), 'kg m^2')
    k = dg.var('k',K_s(I.val,T.av),k_err(T,I), 'N m')
    #beta = dg.var('beta', 2*I.val*gamma.val, beta_err(gamma,T), 'J s')
    #theta_s = StaticG()
    #G_s = dg.var('G_s', G(theta_s.av), G_s_err(k,theta_s,R,M_m,m_m,d), 'N m^2 kg^-2')
    #theta_d = DrivenG()
    #G_d = dg.var('G_d', G(theta_d.val), G_d_err(k,theta_d,R,M_m,m_m,d), 'N m^2 kg^-2')
    
    #dg.printvar(w_b)
    #dg.printvar(I)
    dg.printvar(k)
    #dg.printvar(beta)
    #dg.printvar(theta_s)
    #dg.printvar(G_s)
    #dg.printvar(theta_d)
    #dg.printvar(G_d)
    
    print(G(theta_s.av),G(theta_d.val))
    print(G2(theta_s.av),G2(theta_d.val))
    
    
    return
torun()
