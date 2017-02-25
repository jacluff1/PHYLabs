import numpy as np
import matplotlib.pyplot as plt
import Functions6 as fn

print("working")

# Parameters
CallAll = False # prints and plots all lab results
Test = False

plotA = True
printA = True
saveA = False

figS = fn.figS # figure size for plots
fontS = fn.fontS # plot axis label size
LW = fn.LW
titS = fn.titS # plot title size
binN = fn.binN
N_plot = fn.N_plot
X_diff = fn.X_diff # mili Rad

if CallAll == True:
    plotA = True
    printA = True
    saveA = True
if Test == True:
    plotA = False
    printA = False
    saveA = False
if saveA == True:
    plotA = True
    
# constants
a = fn.a
d  = fn.d
hc = fn.hc

# Error Tracking
e_cnt = fn.e_cnt
dt,e_dt = fn.dt,fn.e_dt # s
red,e_red = fn.red,fn.e_red # nm
green,e_green = fn.green,fn.e_green # nm
L,e_L = fn.L,fn.e_L
e_A,e_V = fn.e_A,fn.e_V
noise_cnt,e_noise_cnt = fn.noise_cnt,fn.e_noise_cnt
noise_V,e_noise_V = fn.noise_V,fn.e_noise_V
V0,e_V0 = fn.V0,fn.e_V0
Cmax,e_Cmax = fn.Cmax,fn.e_Cmax

A_red_params = [0.21317571424375423, 0.4, 198380.42111118653, 146536.85404332707] # mili-volt, mili-rad, nm, nm
B_red_params = [0.37686705355519456, 0.6, 686734.82722076506, 495740.63322611887]
AB_red_params = [1.05, 0.49, 85.0, 67.0, 506]

########################################################################################

def SIX_A(figN):
    C = fn.OpenCSV('01 Stat dark current.csv',1,1)
    X,Y1 = fn.XY('stat',C,0)
    Y2 = fn.CntToVolt(Y1,dt,red,0)
    
    if plotA == True:
        # Dark Current Plot
        CL1 = fn.Leg("data",'b','Cts')            
        CL2 = fn.Leg(0,'k','Cts mean')
        
        VL1 = fn.Leg(0,'b','V(Cts)')
        VL2 = fn.Leg(0,'r','$\delta\ V$')
        VL3 = fn.Leg(0,'k','V(Cts) mean')
        
        fig = plt.figure(figN,figsize=figS)
        
        ax1 = plt.subplot(2,1,1)
        ax1.set_ylabel('CPS',fontsize=fontS)
        ax1.plot(X,Y1,'b-')
        ax1.set_xlim([min(X),max(X)])
        fn.BOX(ax1,X,Y1,'c','k','sig1')
        fn.Legend(ax1,[CL1,CL2],1)
        ax1.set_yticks([np.average(Y1),np.average(Y1)-np.std(Y1),np.average(Y1)+np.std(Y1)])
        ax1.grid(True)
    
        ax2 = plt.subplot(2,1,2)
        ax2.set_xlabel('t (s)',fontsize=fontS)
        ax2.set_ylabel('mV',fontsize=fontS)
        ax2.plot(X,Y2,'b-')
        ax2.set_xlim([min(X),max(X)])
        ax2.plot(X,Y2 + fn.deltaV(Y2,np.std(Y2),dt,e_dt,red,e_red),'r-')
        ax2.plot(X,Y2 - fn.deltaV(Y2,np.std(Y2),dt,e_dt,red,e_red),'r-')
        fn.BOX(ax2,X,Y2,'c','k','sig1')
        fn.Legend(ax2,[VL1,VL2,VL3],1)
        ax2.set_yticks([np.average(Y2),np.average(Y2)-np.std(Y2),np.average(Y2)+np.std(Y2)])
        ax2.grid(True)
    
    # Min/Max Plots
    MIN = fn.OpenCSV('01 Stat minimum counts.csv',1,1)
    MAX = fn.OpenCSV('01 Stat maximum counts.csv',1,1)    
    
    X_MIN,Y_MIN = fn.XY('stat',MIN,-np.average(Y1))
    X_MAX,Y_MAX = fn.XY('stat',MAX,-np.average(Y1))
    V_MIN = fn.CntToVolt(Y_MIN,dt,red,-np.average(Y2))
    V_MAX = fn.CntToVolt(Y_MAX,dt,red,np.average(Y2))
    
    if plotA == True:
        CminL = fn.Leg('data','b','Count Min Data')
        CmaxL = fn.Leg('data','b','Count Max Data')
        VminL = fn.Leg(1,'r','Volt Min')
        VmaxL = fn.Leg(1,'r','Volt Min')
        
        fig2 = plt.figure(figN+1,figsize=figS)
        
        ax1_2 = plt.subplot(2,2,1)
        ax1_2.set_ylabel('Min #/s',fontsize=fontS)
        ax1_2.set_xlabel('t (s)',fontsize=fontS)
        ax1_2.set_xlim([min(X_MIN),max(X_MIN)])
        ax1_2.set_yticks([np.average(Y_MIN),np.average(Y_MIN)-np.std(Y_MIN),np.average(Y_MIN)+np.std(Y_MIN)])
        fn.Legend(ax1_2,[CminL],1)
        ax1_2.scatter(X_MIN,Y_MIN,color='b')
        fn.BOX(ax1_2,X_MIN,Y_MIN,'none','k','sig1')
        
        ax2_2 = plt.subplot(2,2,2)
        ax2_2.set_ylabel('Max #/s',fontsize=fontS)
        ax2_2.set_xlabel('t (s)',fontsize=fontS)
        ax2_2.set_xlim([min(X_MAX),max(X_MAX)])
        ax2_2.set_yticks([np.average(Y_MAX),np.average(Y_MAX)-np.std(Y_MAX),np.average(Y_MAX)+np.std(Y_MAX)])
        fn.Legend(ax2_2,[CmaxL],1)
        ax2_2.scatter(X_MAX,Y_MAX,color='b')
        fn.BOX(ax2_2,X_MAX,Y_MAX,'none','k','sig1')
        
        ax3_2 = plt.subplot(2,2,3)
        ax3_2.set_ylabel('Min mV',fontsize=fontS)
        ax3_2.set_xlabel('t (s)',fontsize=fontS)
        ax3_2.set_xlim([min(X_MIN),max(X_MIN)])
        ax3_2.set_yticks([np.average(V_MIN),np.average(V_MIN)-np.std(V_MIN),np.average(V_MIN)+np.std(V_MIN)])
        fn.Legend(ax3_2,[VminL],1)
        ax3_2.scatter(X_MIN,V_MIN,color='r')
        fn.BOX(ax3_2,X_MAX,V_MIN,'none','k','sig1')
        
        ax4_2 = plt.subplot(2,2,4)
        ax4_2.set_ylabel('Max mV',fontsize=fontS)
        ax4_2.set_xlabel('t (s)',fontsize=fontS)
        ax4_2.set_xlim([min(X_MAX),max(X_MAX)])
        ax4_2.set_yticks([np.average(V_MAX),np.average(V_MAX)-np.std(V_MAX),np.average(V_MAX)+np.std(V_MAX)])
        fn.Legend(ax4_2,[VmaxL],1)
        ax4_2.scatter(X_MAX,V_MAX,color='r')
        fn.BOX(ax4_2,X_MAX,V_MAX,'none','k','sig1')

    if saveA == True:
        fig.savefig('Fig_%s DarkCurrent Stats' % figN)
        fig2.savefig('Fig_%s Max and Min Stats' % (figN+1))
        plt.close(fig)
        plt.close(fig2)
        
        return
    
    CStats,VStats = fn.Stats(Y_MAX),fn.Stats(V_MAX)
    if printA == True:
        print("")
        print("Subtract dark current from results, plot DP for A, B, AB from laser and compare to theory")
        print("The Statistics counter was run with box closed, no bulbs on, no voltage, and outside lights on.")
        print("Figure 1 shows the background noise, or dark current. The top graph shows the Cts. Blue is data, black is mean, cyan represents mean \pm STD")
        print("The cts are converted to volts and error propagtion is used (in red), this didn't really work out so previous standard was used.")
        print("Cnt Stats (noise subtracted): STD = %s, VAR = %s, MEAN = %s" % (CStats[0],CStats[1],CStats[2]))
        print("V Stats (noise subtracted): STD = %s, VAR = %s, MEAN = %s" % (VStats[0],VStats[1],VStats[2]))
    return 

##################################################################################################################################################

def SIX_B(figN):
    
    A = fn.OpenCSV('02 Diff red A.csv',2,1)
    B = fn.OpenCSV('02 Diff red B.csv',1,1)
    AB = fn.OpenCSV('02 Diff red AB.csv',1,1)
        
    XA,YA = fn.XY('diff',A,-noise_V)
    XB,YB = fn.XY('diff',B,-noise_V)
    XAB,YAB = fn.XY('diff',AB,-noise_V)
    
    if plotA == True:

        fig = plt.figure(figN,figsize=figS)
        #fig = plt.figure(figsize=figS)
        
        AL1 = fn.Leg('data','b','Data: A')
        AL2 = fn.Leg(0,'m','Theory: Fit1')
        AL3 = fn.Leg(0,'k','Theory: Fit2')
        ax1 = plt.subplot(3,1,1)
        ax1.set_ylabel('A: mV',fontsize=fontS)
        fn.Legend(ax1,[AL1,AL2,AL3],1)
        ax1.errorbar(XA,YA, xerr=e_A, yerr=e_V, fmt='bs')
        ax1.set_xlim([0,max(X_diff)])
        
        A_red_params = [.2,.4,(a/1000),(red/10)] 
        #A_red_params = [.2,.4,(a/1000),(green/10)]
        ax1.plot(X_diff,fn.SS(X_diff,A_red_params),'m',linewidth = LW/2)
        
        YAFit1,AFitParams1 = fn.SSFit1('A',XA,YA,A_red_params) 
        #ax1.plot(X_diff,YAFit1,'m',linewidth = LW)
        
        YAFit2,AFitParams2 = fn.SSFit2('A',XA,YA,AFitParams1)
        ax1.plot(X_diff,YAFit2,'k',linewidth = LW)   
        
        ABL1 = fn.Leg('data','g','Data: AB')
        ABL2 = fn.Leg(0,'m','Theory: Fit1')
        ABL3 = fn.Leg(0,'k','Theory: Fit2')
        ax2 = plt.subplot(3,1,2)
        ax2.set_ylabel('AB: mV',fontsize=fontS)
        fn.Legend(ax2,[ABL1,ABL2,ABL3],1)
        ax2.errorbar(XAB,YAB, xerr=e_A, yerr=e_V, fmt='gs')
        ax2.set_xlim([0,max(X_diff)])
        
        AB_red_params = [(1.05),.49,(a/1000),(red/10),(d)+100]
        #AB_red_params = [(1.05),.49,(a/1000),(green/10),(d)+10]
        ax2.plot(X_diff,fn.DS(X_diff,AB_red_params),'m',linewidth = LW/2)
        
        YABFit1,ABFitParams1 = fn.DSFit(XAB,YAB,AB_red_params)
        ax2.plot(X_diff,YABFit1,'k',linewidth=LW)
        
        BL1 = fn.Leg('data','r','Data: B')
        BL2 = fn.Leg(0,'k','Theory: Fit1')
        BL3 = fn.Leg(0,'m','Theory: Fit2')
        ax3 = plt.subplot(3,1,3)
        ax3.set_ylabel('B: mV',fontsize=fontS)
        ax3.set_xlabel('$theta [radians]_{-3}$',fontsize=fontS)
        fn.Legend(ax3,[BL1,BL2,BL3],2)
        ax3.errorbar(XB,YB, xerr=fn.e_A, yerr=fn.e_V, fmt='rs')
        ax3.set_xlim([0,max(X_diff)])
        
        B_red_params = [.38,.6,85,60]
        #B_red_params = [.38,.6,(a/1000),(green/10)]
        ax3.plot(X_diff,fn.SS(X_diff,B_red_params),'m',linewidth = LW/2)
        
        YBFit1,BFitParams1 = fn.SSFit1('B',XB,YB,B_red_params) # V0,a,shift,waveL
        #ax3.plot(X_diff,YBFit1,'m', linewidth = LW/2)
        
        YBFit2,BFitParams2 = fn.SSFit2('B',XB,YB,BFitParams1)
        ax3.plot(X_diff,YBFit2,'k', linewidth = LW)
            
    if saveA == True:
        fig.savefig('Fig_%s Diffraction Patterns_green' % (figN))
        plt.close(fig)
    
    if printA == True:
        print("")
        print("Figure %s shows diffraction data for A, B, AB on red laser and fitted using 'SSFit' and 'DSFit'" % figN)
        print("Fit A:")
        fn.PrintParams('A1',AFitParams2)
        print("")
        print("Fit B:")
        fn.PrintParams('B1',BFitParams2)
        print("")
        print("Fit AB:")
        fn.PrintParams('AB',ABFitParams1)
        print("Results:")
        resa,resda = fn.findA(90704.064,92813.117,66827.262)
        print("a +- da = %s +- %s" % (resa,resda))
        
    return

###########################################################################

def SIX_C(figN):
    #C = fn.OpenCSV('09 Stat 10 cts.csv',1,1)
    #X,Y = fn.XY()
    print("single photon diffraction table?")
    return

######################################################################

def SIX_D():
    ans = (1/1000) * (5/1.602E-19) * (670/1239.8)
    ans2 = ans * (1/2.998E8) * (.493)
    print("get 5 mW, convert to eV, divide by energy/photon, gives photons/sec")
    print("photons/sec = (5/1.602E-19) * (1/1000) * (670/1239.8) = %s photons/sec" % ans) 
    print(ans/2.998E8)
    print("multiply photons/sec by (L/c), where L is length of box and c is speed of light.")
    print("this gives number of photons in the box at a time.")
    print("photons in the box at a time = %s" % ans2)
    
    return

############################################################################

def SEVEN_A(figN):
    C = fn.OpenCSV('10 Stat 100 cts.csv',1,1)
    X,Y = fn.XY(0,C,-noise_cnt)
    
    yav = np.average(np.sort(Y))
    
    if printA  == True:
        print("")
        print("Statistis of 100 count Data")
    
    if plotA == True:
        
        fig = plt.figure(figN,figsize = figS)
        ax1 = plt.subplot(2,2,1)
        ax2 = plt.subplot(2,2,2)
        ax3 = plt.subplot(2,2,3)
        ax4 = plt.subplot(2,2,4)
        
        ax1.set_xlim([0,max(X)])
        ax1.set_ylim([0,max(Y)])
        ax1.set_xlabel('t (s)',fontsize=fontS*.75)
        ax1.set_ylabel('$\sigma_1\ Counts$',fontsize=fontS)
        p1_1,p2_1,per_1,hw_1 = fn.BOX(ax1,X,Y,'g','r','sig1')
        ax1.set_yticks([yav, yav-hw_1, yav+hw_1])
        ax1.scatter(X,Y,color='k',s=30)
        
        ax2.set_xlim([0,max(X)])
        ax2.set_ylim([0,max(Y)])
        ax2.set_xlabel('t (s)',fontsize=fontS*.75)
        ax2.set_ylabel('$\sigma_2\ Counts$',fontsize=fontS)
        p1_2,p2_2,per_2,hw_2 = fn.BOX(ax2,X,Y,'c','r','sig2')
        ax2.set_yticks([yav, yav-hw_2, yav+hw_2])
        ax2.scatter(X,Y,color='k',s=30)
        
        ax3.set_xlim([0,max(X)])
        ax3.set_ylim([0,max(Y)])
        ax3.set_xlabel('t (s)',fontsize=fontS*.75)
        ax3.set_ylabel('$68\%\ Counts$',fontsize=fontS)
        p1_3,p2_3,per_3,hw_3 = fn.BOX(ax3,X,Y,'b','r','68')
        ax3.set_yticks([yav, yav-hw_3, yav+hw_3])
        ax3.scatter(X,Y,color='k',s=30)
        
        ax4.set_xlim([0,max(X)])
        ax4.set_ylim([0,max(Y)])
        ax4.set_xlabel('t (s)',fontsize=fontS*.75)
        ax4.set_ylabel('$95\%\ Counts$',fontsize=fontS)
        p1_4,p2_4,per_4,hw_4 = fn.BOX(ax4,X,Y,'m','r','95')
        ax4.set_yticks([yav, yav-hw_4, yav+hw_4])
        ax4.scatter(X,Y,color='k',s=30)
        
        if saveA == True:
            fig.savefig('Fig_%s High Count Statistics' % figN)
            plt.close(fig)
    
    return
    
#############################################################################
    
def SEVEN_B(figN):
    
    C10 = fn.OpenCSV('09 Stat 10 cts.csv',1,1)
    C100 = fn.OpenCSV('10 Stat 100 cts.csv',1,1)
    
    X1,Y1 = fn.XY(0,C10,0)
    X2,Y2 = fn.XY(0,C100,0)
    
    binwidth = 1
    
    if plotA == True:
        
        fig1 = fn.Hist_Plot(figN,X1,Y1,'b',binwidth)
        fig2 = fn.Hist_Plot(figN+1,X2,Y2,'g',binwidth)
        
        fig3 = fn.PLOTS(figN+2,1,'09 Stat 10 cts.csv',False,[np.var(Y1),np.average(Y1)])
        fig4 = fn.PLOTS(figN+3,1,'10 Stat 100 cts.csv',False,[np.var(Y2),np.average(Y2)])
        
        #saveA = False
        if saveA == True:
            fig1.savefig("Fig_%s C10 Scatt Hist" % figN)
            fig2.savefig("Fig_%s C100 Scatt Hist" % (figN+1))
            fig3.savefig("Fig_%s C10 Hist" % (figN+2))
            fig4.savefig("Fig_%s C100 Hist" % (figN+3))
            plt.close(fig1)
            plt.close(fig2)
            plt.close(fig3)
            plt.close(fig4)
            
    return

###############################################################################

def CALL():
    
    if CallAll == True:
        print("")
        SIX_A(1) # figs 1 and 2
        SIX_B(3)
        SIX_C(4)        
        SIX_D()
        SEVEN_A(5)
        SEVEN_B(6)
        
    else:
        print("")
        SIX_A(1)
        SIX_B(3)
        SIX_C(4)
        SIX_D()
        SEVEN_A(5)
        SEVEN_B(6)
        
    return
CALL()


        