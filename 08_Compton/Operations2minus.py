import Functions2minus as fn
import matplotlib.pyplot as plt
import numpy as np

###############################################################################

# CONSTANTS


###############################################################################

def Analysis1(case,plotA,saveA,printA):
    print("Analysis1 running")
    
    if plotA == True:
        
        if case == 1:
            fig = fn.PLOT(case,fn.params0,[fn.calibration],2)
            plt.vlines(fn.binmin, ymin=fn.params0.ymin, ymax=fn.params0.ymax, color='r', linewidth=2)
            plt.vlines(fn.binmax, ymin=fn.params0.ymin, ymax=fn.params0.ymax, color='m', linewidth=2)
            plt.annotate('Ba Channel %s' % fn.binmin, (fn.binmin-10,900), size=20, color='r')
            plt.annotate('Cs Channel %s' % fn.binmax, (fn.binmax-10,1500), size=20, color='m')
            saveas = '01 Calibration bins'
                
        if case == 2:
            fig = fn.PLOT(case,fn.params1,[fn.calibration],2)
            plt.vlines(fn.Ba, ymin=fn.params1.ymin, ymax=fn.params1.ymax, color='r', linewidth=2)
            plt.vlines(fn.Cs, ymin=fn.params1.ymin, ymax=fn.params1.ymax, color='m', linewidth=2)
            plt.annotate('Ba keV %s' % fn.Ba, (fn.Ba-65,900), size=20, color='r')
            plt.annotate('Cs keV %s' % fn.Cs, (fn.Cs-65,1500), size=20, color='m')
            saveas = '02 Calibration keV'
            
        if case == 3:
            fig,AX = fn.PLOT_axis(case,fn.params3,fn.INOUT)
            saveas = '03 With and Without'
            
            
        if case == 4:
            params = fn.params4
            fig = plt.figure(case,figsize=fn.size) 
            saveas = '04 Difference'
            
            for d in fn.IN:
                
                ax = fig.add_subplot(d.axis)
                ax.set_xlabel(params.xlabel,size=fn.FS)
                ax.set_ylabel(params.ylabel,size=fn.FS)
                ax.set_xlim([params.xmin,params.xmax])
                ax.set_ylim([params.ymin,params.ymax])
                ax.set_title(d.label[:6])
                
                data = fn.OpenCSV(d.CSV,params.LL,int)
                X,Y = fn.XY(data,params.convert)
                
                d2 = fn.OUT[d.index]               
                data2 = fn.OpenCSV(d2.CSV,params.LL,int)
                X2,Y2 = fn.XY(data2,params.convert)
                
                i_max1,i_max2 = len(Y),len(Y2)
                if i_max1<= i_max2:
                    i = i_max1
                else:
                    i = i_max2
                
                Ydiff = Y[:i]-Y2[:i]
                ypeak = max(Ydiff)
                xpeak = X[(np.abs(Ydiff-ypeak)).argmin()]
                            
                plt.plot(X[:i],Ydiff, label='Difference')
                ax.vlines(xpeak, ymin=params.ymin, ymax=ypeak, color='r', linewidth=2)
                ax.hlines(ypeak, xmin=params.xmin, xmax=xpeak, color='r', linewidth=2)
                ax.vlines(fn.Cs, ymin=params.ymin, ymax=ypeak, color='m', linewidth=2)
                ax.hlines(ypeak, xmin=xpeak, xmax=fn.Cs, color='m', linewidth=2)
                ax.annotate('( %s , %s )' % (xpeak,ypeak), (params.xmin+50,ypeak+20), fontsize=20, color='r')
                ax.annotate('$\Delta\ E$ = %s keV' % (fn.Cs-xpeak), (params.xmax-350,ypeak+20), fontsize=20, color='m')
                
                rewrite = True
                if rewrite == True:
                    np.save(fn.DIFF[d.index].NP+'Y.npy', Ydiff)
                    np.save(fn.DIFF[d.index].NP+'X.npy', X[:i])
                    
            plt.tight_layout()
        
        if case == 5:
            params = fn.params5  
            saveas = '05 ComptonFit'
            
            plt.close()
            plt.close()
            fig = plt.figure(case,figsize=fn.size)
            plt.xlabel(params.xlabel,size=fn.FS)
            plt.ylabel(params.ylabel,size=fn.FS)
            plt.xlim([params.xmin,params.xmax])
            plt.ylim([params.ymin,params.ymax])
            
            X,Y = [],[]
            X_err,Y_err = [],[]
            UNC = []
            for d in fn.DIFF:
                theta = np.deg2rad(d.deg)
                x = 1-np.cos(theta)
                y = 1/d.fit_alg[1]
                x_err = np.sin(theta)*.1
                y_err = 7.262/d.fit_alg[1]**2
                X.append(x)
                Y.append(y)
                X_err.append(x_err)
                Y_err.append(y_err)
                unc = np.sqrt( (x_err*(y-(1/661))**(-1))**2 + ( (x*y_err*661**2)/((661*y - 1)**2) )**2 )
                UNC.append(unc)
            uncertainty = np.average(UNC)
                
            X2,Y2,m,E = fn.ComptonFit(X,Y)
            Y3 = fn.Eq2(X2,511,fn.Cs)
            
            #data
            #plt.scatter(X,Y,c='b',marker='s',s=50)
            plt.errorbar(X,Y, xerr=X_err, yerr=Y_err, c='b', fmt='bs', markersize=10, label='data')
            #fit
            plt.plot(X2,Y2,'r',linewidth=2, label='fit')
            plt.annotate('FIT: y = (x/%s) + (1/%s)' % (round(m,2),round(E,2)), (.18,.0026), fontsize=20, color='r')
            #benchmark
            plt.plot(X2,Y3,'m',linewidth=2, label='benchmark')
            plt.annotate('Benchmark: y = (x/%s) + (1/%s)' % (511,fn.Cs), (.4,.0021), fontsize=20, color='m')
            plt.legend(loc=2)
                
                         
        if saveA == True:
            fig.savefig(fn.Gpath + 'Analysis1\\' + saveas + '.png')
            plt.close()
        
    if printA == True:
        print("")
        print("Analysis Part 1")
        print("Fitting Channels = %s , %s" % (fn.binmin,fn.binmax))
        print("Fitting keV = %s , %s" % (fn.Ba,fn.Cs))
        print("A , B = %s , %s" % (fn.A,fn.B))
        print("Compton Fit: m_e c^2 , E = %s , %s" % (m,E))
        err1,err2 = (511-m),(fn.Cs-E)
        print("Compton Error: %s , %s" % (err1,err2))
        print("Compton percent Error: %s , %s" % ( (err1*100/511),(err2*100/fn.Cs) ))
        print("Uncertainty in electron rest mass = %s" % uncertainty)
        print("")
    
    print("Analysis1 complete")
    return

def Analysis2(case,plotA,saveA,printA):
    print("Analysis2 running")

    if case == 1:
        params = fn.params4
        plt.close()
        plt.close()
        fig = plt.figure(case,figsize=fn.size)
        saveas = '01 GaussFit'
        
        X_ptr,Y_ptr = [],[] # peak to total ratios
        for d in fn.DIFF:
            
            ax = fig.add_subplot(d.inA.axis)
            ax.set_xlabel(params.xlabel,size=fn.FS)
            ax.set_ylabel(params.ylabel,size=fn.FS)
            ax.set_xlim([0,params.xmax])
            ax.set_ylim([params.ymin,params.ymax])
            ax.set_title(d.inA.label[:6])
            
            X = np.load(d.NP + 'X.npy')
            Y = np.load(d.NP + 'Y.npy')
            Y = (Y/d.QE)-min(Y)
            LL,UL,fit_params = fn.GaussFit(X,Y,d)
            Y_g = fn.gaus(X,*fit_params)
            
            ax.plot(X,Y,'b', linewidth=.5, label='data')
            ax.plot(X,fn.gaus(X,*d.fit_manual),'r',linewidth=1.5, label='manual fit')
            ax.plot(X,Y_g,color='k',linewidth=2, label='algorithm fit')
            
            for i in range(len(X)):
                x,y = X[i],Y_g[i]
                if y>d.fit_alg[0]*.001:
                    width=d.fit_alg[1]-x
                    break
            LLg,ULg = d.fit_alg[1]-width,d.fit_alg[1]+width # energy limits of gaussian region
            ax.vlines([LLg,ULg], ymin=0, ymax=params.ymax, color='m',linewidth=1.5,label='gauss region')
            #print(d.index,width)
            ax.legend(loc=2,fontsize=20)
            
            rewrite = True
            if rewrite == True:
                LLti = (np.abs(X-100)).argmin()
                LLgi,ULgi = (np.abs(X-LLg)).argmin(), (np.abs(X-ULg)).argmin()
                Xgr,Ygr = X[LLgi:ULgi],Y_g[LLgi:ULgi]
                Xtr,Ytr = X[LLti:ULgi], Y[LLti:ULgi]
                
                peak,total = 0,0
                for i in range(len(Xgr)):
                    peak += Ygr[i]
                for i in range(len(Xtr)):
                    total += Ytr[i]
                Y_ptr.append(peak/total)
                X_ptr.append(d.peakx)
                np.save(fn.NPpath + 'ptr_X.npy',X_ptr)
                np.save(fn.NPpath + 'ptr_Y.npy',Y_ptr)
                print(peak/total)

        plt.tight_layout()
                
    if case == 2:
        saveas = '02 Differential Cross Section'
        X_yield,Y_yield = [],[]
        err = []
        X = np.linspace(0,100,1000)
        for d in fn.DIFF:
            Y = np.load(d.NP+'Y.npy')
            Y = Y * (d.PtT/d.QE)
            X_yield.append(d.deg)
            y = np.sum(Y)*1e27
            d.CS = fn.diff_crossSection(y)
            Y_yield.append(d.CS)
            err.append(fn.diff_crossSection(y,ret='two')/150)
        err = np.array(err).astype(float)
        Err = np.average(err)
        
        X_mel = [20,30,40,60,80,100]
        Y_mel = [55.2,42.9,35.8,23.9,18.0,15.8]
            
        #plt.close()
        fig = plt.figure(case,fn.size)
        plt.xlabel('Scattering Angle Degrees',fontsize=20)
        plt.ylabel('$d\\sigma/d\\Omega\ [10^{-27}\ cm^2/sr]$',fontsize=20)
        plt.xlim([0,100])
        #plt.scatter(X_yield,Y_yield, c='b', s=30, marker='s', label='data')
        plt.errorbar(X_yield,Y_yield, xerr=.5, yerr=err, color='b', markersize=8, fmt='bs', label='data')
        plt.plot(X,fn.KleinNishina(X), color='r', label = 'Klein-Nishina')
        plt.scatter(X_mel,Y_mel, color='m', s=30, marker='s', label="Melissianos")
        plt.legend(loc=3,fontsize=20)
        plt.annotate("Its practically the same line!!", (50,40), color='r',fontsize=20)
         
    if saveA == True:
        fig.savefig(fn.Gpath + 'Analysis2\\' + saveas + '.png')
        plt.close()
    
    if printA == True:
        print("")
        print("Analysis2")
        print("Manual fit params: [ height , average , sigma]")
        for d in fn.DIFF:
            print(d.fit_manual)
        print("Alg fit params: [ height , average , sigma]")
        for d in fn.DIFF:
            print(d.fit_alg)
        print("Differential Cross Sections: 1e-27 cm^2")
        for d in fn.DIFF:
            print(d.CS)
        print("error = %s 1e-27 cm^2" % Err)
        print("")
    
    print("Analysis2 completed")
    return

            
            

###############################################################################
