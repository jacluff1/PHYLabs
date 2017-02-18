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
        
# GLOBAL VARIABLES-------------------------------------------------------------

nppath = "../npy/"
csvp = "../csv/"
gp = "../graphs/"
per = '$\%$'

csvs = ["angle_cal","static1","static2","driven","equilib","extra_decay"]
npys = csvs

# MEASUREMENTS-----------------------------------------------------------------

# m and M
M1 = ball(1054.5,5.15,True)
M2 = ball(1052.9,5.15,True)
m1 = ball(14.5,.85,True)
m2 = ball(14.5,.85,True)
M = ball((M1.m+M2.m)/2,(M1.d+M2.d)/2,False)
m = ball((m1.m+m2.m)/2,(m1.d+m2.d)/2,False)

# angle calibration
S = dg.var('S',np.array([43.6,28.4,18.3,11.9,7.8,5.1])/100, .5/100, 'meters')
L = dg.var('L',762/100, .1, 'meters')

# boom
l_b = dg.var('l_b',14.9/100, .1/100, 'm') # length of boom
m_b = dg.var('m_b',7.1/1000, .05/1000, 'kg') # mass of boom
w_b = dg.var('w_b',.8/100, .05/100, 'm') # width of boom

# cavendish box
w_t = dg.var('w_t',.1/100, .05/100, 'm') # glass thickness
W = dg.var('W',2.5/100, .05/100, 'm') # width of box

# m and m'
m_m = dg.var('m_m',m.m, .05/1000, 'kg') # mass of m
m_d = dg.var('m_d',m.d, .05/100, 'm') # diameter of m
m_r = dg.var('m_r',m_d.val/2, m_d.err/2, 'm') # radius of m

# M and M'
M_m = dg.var('M_m',M.m, .05/1000, 'kg') # mass of M
M_d = dg.var('M_d',M.d, .05/100, 'm') # diameter of M
M_r = dg.var('M_r',M_d.val/2, M_d.err/2, 'm') # radius of M

# FITTING PARAMATERS-----------------------------------------------------------

par_m = [0,10,1.5,(2*np.pi)/(.385)]
theta_e = dg.var('theta_e',0.0, 0.1, 'mRad') # equilibrium angle
A = dg.var('A',10.5, 0.2, 'mRad') # initial amplitude of decaying wave
gamma = dg.var('gamma',1.6/1000, 0.1/1000, 'rad s^-1')
omega = dg.var('omega',16.3/1000, 0.2/1000, 'rad s^-1') # angular frequency 

T = dg.var('T',[396., 385., 380., 380., 382., 382.5], 'std', 's') # Period
C = dg.var('C',.89,.01,None)

# LAB EQUATIONS ---------------------------------------------------------------

def d(l_b,m_d): # axis to small spheres - meters
    val = (1/2) * (l_b.val - m_d.val)
    err = (1/2) * np.sqrt( l_b.err**2 + m_d.err**2 )
    return dg.var('d',val,err,'m')

def R(W,w_t,M_d): # distance between m and M - meters
    val = (1/2) * (W.val + M_d.val) + w_t.val
    err = (1/2) * np.sqrt( W.err**2 + M_d.err**2 + 4*w_t.err**2 )
    return dg.var('R',val,err,'m')

def I(m_m,d,m_r,m_b,l_b,w_b): # moment of inertia - kg m^2
    val1 = (2*m_m.val) * (d.val**2 + (2/5)*m_r.val**2)
    val2 = (m_b.val/12) * (l_b.val**2 + w_b.val**2)
    val = val1 + val2
    
    var1 = (2*d.val**2*m_b.err)**2 + (4*m_m.val*d.val*d.err)**2 + ((4/5)*m_r.val**2*m_m.err)**2 + ((8/5)*m_m.val*m_r.val*m_r.err)**2
    var2 = ((1/12)*l_b.val**2*m_b.err)**2 + ((1/6)*m_b.val*l_b.val*l_b.err)**2 + ((1/12)*w_b.val**2*m_b.err)**2 + ((1/6)*w_b.val*m_b.val*w_b.err)**2
    err = np.sqrt(var1 + var2) 
    return dg.var('I',val,err,'kg m^2')

def K(omega,gamma,I): # torsion constant - N m
    val = I.val * (omega.val**2 + gamma.val**2)
    
    var1 = ((omega.val**2 + gamma.val**2)*I.err)**2
    var2 = (2*I.val*omega.val*omega.err)**2
    var3 = (2*I.val*gamma.val*gamma.err)**2
    err = np.sqrt( var1 + var2 + var3 )
    return dg.var('k',val,err,'N m')

def G(theta,k,R,M_m,m_m,d,static=True): # gravitational constant - N m^2 kg^-2
    if static == True:
        th_err = theta.err/1000 # mRad -> rad
        th_val = theta.av/1000 # mRad -> rad
    else:
        th_err = theta.err/1000
        th_val = theta.val/1000
    
    m_h = .34/1000 # mass missing from boom where small spheres are - kg
    f_d = 3.5/100
    f_b = 0.19
    m_factor = (m_m.val - m_h)*(1-f_d) + (m_b.val*f_b)
    
    #val = (k.val * th_val * R.val**2)/(2 * M_m.val * m_m.val * d.val)
    val = (k.val * th_val * R.val**2)/(2 * M_m.val * m_factor * d.val) # corrected
    
    var1 = (k.err/k.val)**2
    var2 = (th_err/th_val)**2 
    var3 = (2*R.err/R.val)**2 
    var4 = (m_m.val*m_m.err)**2 
    var5 = (M_m.val*M_m.err)**2 
    var6 = (d.val*d.err)**2
    err = val * np.sqrt( var1 + var2 + var3 + var4 + var5 + var6 )
    return dg.var('G',val,err,'N m^2 kg^-2') 

def laser(theta_comp,C,theta_comp_err=.05): # theta_laser - mRad
    if type(theta_comp) != float:
        theta_comp = np.array(theta_comp)
        theta_comp = np.average(theta_comp)
    val = theta_comp*C.val
    err = np.sqrt( (theta_comp_err/theta_comp)**2 + (C.err/C.val)**2 )
    return dg.var('laser',val,err,'mRad')

def theta_laser_err(theta,theta_err,C):
    return np.sqrt( (C.val*theta_err/theta)**2 + (theta*C.err/C.val)**2 )

# DATA-------------------------------------------------------------------------

def Write(): # writes all .csv files to .npy files
    for f in csvs:
        dg.CSV(csvp,nppath,f)

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

def theta_driven(X,Y,dx,N=13,printA=True):
    
    tmaxima = np.array([.53, .92, 1.31, 1.70, 2.08, 2.50])
    tminima = np.array([.31, .73, 1.12, 1.51, 1.89, 2.28, 2.67])
    extrema = dg.extrema(X,Y,dx,tminima,tmaxima,ret='tog')
    theta = extrema[:,1]
    laser_err = theta_laser_err(np.average(theta),.05,C)
    if printA == True:
        print("")
        print(theta)
        print(len(theta))
    
    def thetasum(ni,nf):
        LL,UL = ni-1, nf
        total,n = 0,0
        for t in theta[LL:UL]:
            theta_n = t * (-1)**n 
            total += theta_n
            n += 1
        return total
    
    def X1(N,printA=False): # equation 8
        val1 = theta[1-1] - theta[N-1]
        val2 = thetasum(1,N-1)
        val = 1 - val1/val2
        
        err1 = laser_err * (1-val)
        err2 = np.sqrt( (N-1)*(1-val)**2 + (2*val) )
        err3 = abs(theta[1-1] - theta[N-1])
        err = err1 * err2/err3
        
        if printA == True:
            print("x1 = %s +/- %s" % (val,err) )
        
        return dg.var('x1',val,err,None)
    
    def X2(N,printA=False):# equation 8a
        val1 = theta[2-1]-theta[N-1-1]
        val2 = thetasum(2,N-2)
        val = 1 - val1/val2
        
        err1 = 1/abs(theta[1-1]-theta[N-1])
        err2 = np.sqrt( (N-1)*(1-val)**2 + (2*val) )
        err = laser_err * (1-val) * err1 * err2

        if printA == True:
            print("")
            print("x2 = %s +- %s" % (val,err))
        return dg.var('x2',val,abs(err),None)
    
    def X3(N,printA=False):#equation 7
        Xs = []
        for t in theta[:-2]:
            i = dg.nearest(theta,t)
            x = -(theta[i+2] - theta[i+1])/(theta[i+1] - theta[i])
            Xs.append(x)
        Xs = np.array(Xs)
        x = np.average(Xs)
        
        if printA == True:
            print("")
            print("x3 = %s" % Xs)
            print("x3 = %s +- %s" % (x,np.std(Xs)))
        return dg.var('x3',x,np.std(Xs),None)
    
    def X4(gamma,T,printA=False):
        val = np.exp(-gamma.val*T.val/2)
        val = np.average(val)
        err = (val/2) * np.sqrt( (T.val*gamma.err)**2 + (gamma.val*T.err)**2 )
        err = np.average(err)
        if printA == True:
            print("x4 = %s +- %s" % (val,err))
        return dg.var('x',val,err,None)
    
    def thetaD(N,printA=False): # equaiton 15
        val1 = (1-x.val)
        val2 = thetasum(1,N)
        val3 = -theta[1-1] + (x.val*theta[N-1])
        val4 = (N-1)*(1+x.val)
        val = (val1*val2 + val3)/val4
        
        err1_1 = theta_laser_err(np.average(theta),.05,C)
        err1_2 = np.sqrt( (N-1)*(1-x.val)**2 + 2*x.val )
        err1 = err1_1 * err1_2/( (N-1)*(1+x.val) )
        
        err2_1 = x.err/( (N-1)*(1+x.val)**2 )
        err2_2 = 2*thetasum(1,N-1) + theta[N-1] - theta[1-1]
        err2 = err2_1 * err2_2
        
        err = np.sqrt( err1**2 + err2**2 )
        
        if printA == True:
            print("")
            print("thetad val1 = %s" % val1)
            print("thetad val2 = %s" % val2)
            print("thetad val3 = %s" % val3)
            print("thetad val4 = %s" % val4)
            print("thetad val = %s" % val)
        
        return dg.var('theta_d',val,err,'mRad')
    
    def thetaD2(N,printA=False): # equation 14
        displ = []
        for t in theta[:-2]:
            i = dg.nearest(theta,t)
            d1 = (-1)**(i+1)/( (1+x.val) )
            d2 = (x.val * theta[i]) + (1-x.val)*theta[i+1] - theta[i+2]
            displ.append(d1 * d2)
        displ = np.array(displ)
        av = np.average(displ)
        #err = np.std(displ)
        err1_1 = theta_laser_err(np.average(theta),.05,C)
        err1_2 = np.sqrt( (N-1)*(1-x.val)**2 + 2*x.val )
        err1 = err1_1 * err1_2/( (N-1)*(1+x.val) )
        
        err2_1 = x.err/( (N-1)*(1+x.val)**2 )
        err2_2 = 2*thetasum(1,N-1) + theta[N-1] - theta[1-1]
        err2 = err2_1 * err2_2
        
        err = np.sqrt( err1**2 + err2**2 )
        
        if printA == True:
            print("")
            print("thetad 2 = %s" % displ)
            print("av = %s" % av)
        return dg.var('thetaD 2',av,err,'mRad')
    
    def ThetaD3():
        Q = omega.val/(2*gamma.val)
        def the(t):
            return theta_e.val + (4*Q/np.pi)*theta_s.val*(1-np.exp(-gamma.val*t))*np.cos(omega.val*t)
        thetas = the(extrema[:,0])
        return dg.var('thetaD 3',thetas,'std','mRad')
    N=13
    if N == 11:
        theta = theta[:-2]  
        #x = X1(11,printA=True)
        #x = X2(11,printA=True)
        #x = X3(11,printA=True)
        x = X4(gamma,T)
        #x = dg.var('x_test',1.4,.1,None)
        theta_d = thetaD2(11) 
    else:
        #x = X1(13,printA=True)
        #x = X2(13,printA=True)
        #x = X3(13,printA=True)
        #x = X4(gamma,T,printA=True)
        x = dg.var('x_test',.74,.01,None)
        #theta_d = thetaD2(13,printA=True)
        theta_d = ThetaD3()
        
    return extrema, theta_d
    
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

def DrivenG(N=13,plot=False,printA=False):
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
    
    ext,theta_d = theta_driven(X,Y,6,N=N)
    
    d2_ext = dp.data(ext[:,0],ext[:,1],'*','g','$\\theta_n$')
    
    if plot==True:
        ax = dp.ax([data,d1,d2_env,d2_man,d2_fit,d2_ext],111,'time [$10^3$ sec]','$\\theta_{laser}$ [mRad]','Driven Measurement')
        ax.notes = [note]

        dp.plot([ax], name=gp+'0_4_DrivenG')
    
    if printA == True:
        print("fitting parameters")
        print(par_f)
    
    return theta_d
#DrivenG(printA=True,plot=True)
    
d_var = d(l_b,m_d)
R_var = R(W,w_t,M_d)
I_var = I(m_m,d_var,m_r,m_b,l_b,w_b)
k_var = K(omega,gamma,I_var)  
theta_s = StaticG()
theta_d = DrivenG(N=11)
G_s = G(theta_s,k_var,R_var,M_m,m_m,d_var)
G_d = G(theta_d,k_var,R_var,M_m,m_m,d_var,static=False)

def Results():

    print("")
    dg.printvar(d_var)
    dg.printvar(R_var)
    dg.printvar(I_var)
    dg.printvar(k_var)
    print("")
    dg.printvar(theta_s)
    dg.printvar(theta_d)
    print("")
    dg.printvar(G_s)
    dg.printvar(G_d)
    print("")
Results()
print(G_d.val)

    
  