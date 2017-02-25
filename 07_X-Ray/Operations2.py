import Functions1 as fn1
import Functions2 as fn
import numpy as np
import matplotlib.pyplot as plt

###############################################################################

# Parameters

testA = True

###############################################################################

def Part(case,plotA,saveA,printA):
    
    # plot energy spectrum, identify the peaks
    # use theory equations to plot a "vline" at the degrees where these peaks would occur and see how they line up
    # Compute K_alpha and K_beta and compare to known values.
    # Compute K_{2,3} (and compare?)
    
    plt.close()
    plt.close()
    fig1 = plt.figure(case,figsize=fn1.size)
    fig2 = plt.figure(case+1,figsize=fn1.size)
    
    n  = np.array([1,2,3])
        
    theta_alpha = fn.Bragg(n,fn1.Mo.alpha)
    theta_beta = fn.Bragg(n,fn1.Mo.beta)
    
    bohr23 = fn1.Eq1(74,2,3)
    theta23 = fn.Bragg(n,bohr23)
    
    if plotA == False:
        plt.close(fig1)
        plt.close(fig2)
        
    for p in fn.parts:
        
        LL,UL = 9,285
        
        data = fn.OpenCSV(fn.path+p.CSV,LL,UL)
        Xd,Xr,Y = fn.XY(data,0,1)
        
        theta_U = fn.Bragg(1,p.U)
        
        if plotA == True:
        
            xticks = np.arange(2.5,30+2.5,2.5)
            ax = fig1.add_subplot(p.plot[0])
            ax.set_xlabel('degrees',fontsize=fn1.FS)
            ax.set_ylabel('CPS',fontsize=fn1.FS)
            ax.set_xlim([min(Xd),max(Xd)])
            ax.set_xticks(xticks)
            
            dL = fn1.Leg(0,p.plot[1],p.legend[1])
            UL = fn1.Leg(0,'k','U = 35')   
            
            if p.source == 'Mo':
                
                theta1 = theta_alpha
                theta2 = theta_beta
                ax.vlines(theta2, ymin=0, ymax=max(Y), color='c', linewidth=2)
                aL = fn1.Leg(0,'m','$K_\\alpha$ lines')
                bL = fn1.Leg(0,'c','$K_\\beta$ lines') 
                fn1.Legend(ax,[dL,aL,bL,UL],p.legend[0])
                ax.plot(Xd-.75,Y,p.plot[1])
                
            else:
                
                theta1 = theta23
                theta2 = 'none'
                aL = fn1.Leg(0,'m','$K_{2,3}$')
                fn1.Legend(ax,[dL,aL,UL],p.legend[0])
                ax.plot(Xd,Y,p.plot[1])
                
            ax.vlines(theta_U, ymin=0, ymax=max(Y), color='k', linewidth=3)
            ax.vlines(theta1, ymin=0, ymax=max(Y), color='m', linewidth=2)
            
    if printA == True:
        print("")
        print("[ Mo ]")
        print("theta_alpha = %s" % theta_alpha)
        print("theta_beta =  %s" % theta_beta)
        print("")
        print("[ W74 ]")
        print("theta_23 = %s" % theta23)
    
    #figure 2
    data1 = fn.OpenCSV(fn.path+fn.II1.CSV,14,110)
    X1d,X1r,Y1 = fn.XY(data1,0,1) 
    
    data2 = fn.OpenCSV(fn.path+fn.II2.CSV,9,105)
    X2d,X1r,Y2 = fn.XY(data2,0,1)
    
    xticks = np.arange(3,12.5+2.5,2.5)
    ax = fig2.add_subplot(111)
    ax.set_xlabel('degrees',fontsize=fn1.FS)
    ax.set_ylabel('CPS',fontsize=fn1.FS)
    ax.set_xlim([3,12.5])
    ax.set_xticks(xticks) 
    
    ax.plot(X1d,Y1,fn.II1.plot[1], X2d,Y2,fn.II2.plot[1])

    L1 = fn1.Leg(0,fn.II1.plot[1],'U = 35')
    L2 = fn1.Leg(0,fn.II2.plot[1],'U = 30')
    fn1.Legend(ax,[L1,L2],1)
        
    if saveA == True:
        fig1.savefig("Fig_6.1 Bragg.png")
        fig2.savefig("Fig_6.2 Monochromator.png")
        plt.close(fig1)
        plt.close(fig2)
        
    return
        
        