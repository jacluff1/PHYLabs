import Functions1 as fn
import numpy as np
import matplotlib.pyplot as plt

##################################################################################################

# PARAMETERS

testA = False

##################################################################################################

def PartIII(case,plotA,saveA,printA): # there are 3 cases
    if testA == True:
        print("PartIII Running")
        plt.close
        
    # part 5: record spectra for all elements in carousel
    
    # copy path extension to the "Carousel" folder and add an extra slash between subfolders
    path = "C:\\Users\\Jacob/Google Drive\\School\\Physics Degree\\Current Classes\\PHY 334 Advanced Lat 1\\X-Ray\\Excel Files\\Carousel\\"     
    
    # add legends to the "Carousel" class systems' parameters
    for el in fn.Samples:
        leg = fn.Leg('data',el.params[1],el.name)
        el.params[3] = leg
    
    if plotA == True:
        
        if case == 1:
            table = fn.OpenCSV(path+fn.Cu.title,4,2052,int) # full range
            X,Y = fn.XY(False,table,0,1) # X,Y = channel,counts
            Cu_alpha,Am241 = 270,1970
            xticks = [Cu_alpha,Am241]
                
            fig = plt.figure(case,figsize=fn.size)
            plt.xlabel('Channel',fontsize=fn.FS)
            plt.ylabel('Counts',fontsize=fn.FS)
            plt.xticks(xticks)
            plt.vlines([Cu_alpha,Am241,268,272,1968,1972], ymin = min(Y), ymax = max(Y))
            plt.plot(X,Y,'r')
            plt.title('Copper Spectrum: Finding A,B',fontsize=(fn.FS+4))
            
            if printA == True:
                print("")
                print("Part III , Case 1:")
                print("Compare Cu_alpha and 124 Am lines to find A,B parameters.")
                print("Fitting channels = %s , %s" % (Cu_alpha,Am241))
                print("Adjusting paramaters of figure 1 manually, error in fitting channels = \pm %s" % 1)
                print("Use 'keV = A + Bx', the fitting channels, and Cu_alpha and 241Am to find A and B")
                print("A , B = %s , %s" % (fn.A,fn.B))
                print("Using last significant figure calculated : dA , dB = %s , %s" % (fn.dA,fn.dB))
                print("using error propagations, error in keV (x-axis) is found using fn.d_keV(keV)")
                print("")
    
    if case != 1:
        
        fig = plt.figure(case,figsize=fn.size)          
        for el in fn.Samples:
            
            # part 6: fit peak positions to Mosley's law (eq 2).
            # explain the physical meaning of the fit values
            # part 7: Compare measured energy for the Cu K_alpha line to Bohr model value.
                
            if case == 2:
                table = fn.OpenCSV(path+el.title,4,2052,int) # full range
                X,Y = fn.XY(True,table,0,1) # X,Y = keV,counts
                fn.plotAxis(fig,el.params[0],'plot',X,Y,el,[el.alpha,59.5])
                
            if case == 3:
                LL,UL = fn.KeVToBin(fn.A,fn.B,el.ROI[0]),fn.KeVToBin(fn.A,fn.B,el.ROI[1])
                table = fn.OpenCSV(path+el.title,LL+4,UL+4,int) # ROI only
                X,Y = fn.XY(True,table,0,1) # X,Y = keV,counts
                
                kev = fn.Eq2(el.Z,el.C,el.sigma)
                fn.plotAxis(fig,el.params[0],'plot',X,Y,el,[kev,el.beta,el.alpha])
                plt.vlines(kev, ymin=min(Y), ymax=max(Y), color='orange',linewidth=2)
                
                MosleyL = fn.Leg(0,'orange','Mosley Fit')
                fn.Legend(plt,[el.params[3],MosleyL],el.loc)
                
                if printA == True:
                    print("")
                    print("%s:" % el.name)
                    print("Plotted range : %s , %s keV" % (el.ROI[0],el.ROI[1]))
                    print("Mosley Fit: C , sigma , keV , err_keV = %s , %s , %s , %s" % (el.C,el.sigma,kev,el.err_mosley))
            
            if testA == True: break
                
    if printA == True:
        
        if case == 2:
            print("")
            print("Part III , Case 2:")
            print("Fig 3.2 shows the full spectrum of each element in the wheel")
            print("the alpha line for the element and the  241 Am line is plotted to ensure that they line up.")
            print("")
        
        if case == 3:
            print("")
            print("Part III , Case 3:")
            print("only region of interests are plotted")
            print("Mosley Fit C : %s +/- %s" % (np.average(fn.mos_C),np.std(fn.mos_C)))
            print("Mosley Fit Sigma : %s +/- %s" % (np.average(fn.mos_Sigma),np.std(fn.mos_Sigma)))
            print("Mosley Fit Error : %s" % (np.average(fn.mos_error) + np.std(fn.mos_error)))
            print("Need to explain physcal meaning of Mosley fit values")
            Cu_boar = fn.Eq1(fn.Cu.Z,1,2)
            print("Boar Energy of Cu: %s keV" % Cu_boar)
            print("Cu K_alpha : %s keV" % fn.Cu.alpha) 
            print("Boar:K_alpha = %s" % (Cu_boar/fn.Cu.alpha))
            print("")
        
    if saveA == True:
        if case == 1:
            fig.savefig('Fig_3.1 Find A,B.png')
        if case == 2:
            fig.savefig('Fig_3.2 Full Spectrum.png')
        if case == 3:
            fig.savefig('Fig_3.3 Mosley\'s Law.png')
        plt.close(fig)
    
    if testA == True:
        print("PartII completed")
    return

def PartIV(case,plotA,saveA,printA): # Spin-Orbit K_beta doublet
    if testA == True:
        print("PartIV started")
        plt.close()
    
    # 8 find position and intensities, constrain widths of spin-orbit pair to be equal
    # 9 find the power law for spin-orbit splitting (Eq3)
    
    # copy path extension to the "Carousel" folder and add an extra slash between subfolders
    path = "C:\\Users\\Jacob/Google Drive\\School\\Physics Degree\\Current Classes\\PHY 334 Advanced Lat 1\\X-Ray\\Excel Files\\Carousel\\"     
    
    if plotA == True:
        
        fig = plt.figure(case,figsize=fn.size)
        
        for el in fn.Samples:
            dataL = fn.Leg(0,el.params[1],el.name)
            doubL = fn.Leg(0,'orange','Doublet')
            
            if case == 1:
                LL,UL = fn.KeVToBin(fn.A,fn.B,el.SO_ROI[0]),fn.KeVToBin(fn.A,fn.B,el.SO_ROI[1])
                table = fn.OpenCSV(path+el.title,LL+4,UL+4,int) # ROI only
                X,Y = fn.XY(True,table,0,1) # X,Y = keV,counts
                fn.plotAxis(fig,el.params[0],'plot',X,Y,el,[el.beta,min(X),max(X)])
                
                if el.doublet == 'no':
                    fn.Legend(plt,[dataL],1)
                
                if el.doublet != 'no':
                    x,y = el.doublet
                    plt.vlines(x, ymin=0, ymax=y, color='orange', linewidth=2)
                    plt.hlines(y, xmin=0, xmax=x ,color='orange', linewidth=2)
                    fn.Legend(plt,[dataL,doubL],1)
                
                if printA == True:
                    print("%s : [Doublet pos,cps] , [Doublet err] = %s , %s" % (el.name,el.doublet,el.err_doub))
            
            if case == 2:
                if el.doublet != 'no':
                    LL,UL = fn.KeVToBin(fn.A,fn.B,el.SO_ROI[0]),fn.KeVToBin(fn.A,fn.B,el.SO_ROI[1])
                    table = fn.OpenCSV(path+el.title,LL+4,UL+4,int) # ROI only
                    X,Y = fn.XY(True,table,0,1) # X,Y = keV,counts
                    
                    ax = fig.add_subplot(el.doub_axis)
                    ax.plot(X,Y,el.params[1])
                    ax.set_xlim([min(X),max(X)])
                    ax.set_xlabel(fn.Xlabel,fontsize=fn.FS)
                    ax.set_ylabel(fn.Ylabel,fontsize=fn.FS)
                    
                    n = el.doubFit[0]
                    doub = fn.Eq3(el.Z,n)
                    ax.vlines(doub+el.beta, ymin=0, ymax=max(Y), color='orange', linewidth=2)
                    ax.vlines([el.D_ROI[0],el.D_ROI[1]], ymin=0, ymax=max(Y), color='r', linewidth=2)
                    ax.vlines(el.beta, ymin=0, ymax=max(Y),color='g', linewidth=2)                    
                    xticks = [el.beta,el.D_ROI[0],el.D_ROI[1]]
                    ax.set_xticks(xticks)
                    
                    beta = fn.Leg(0,'g','$K_\\beta$')
                    fitL = fn.Leg(0,'orange','$Z^{%s}$' % n)
                    fn.Legend(plt,[dataL,beta,fitL],1)
                    
                    if printA == True:
                        print("n , E_doublet , err_E_doublet_per = %s , %s , %s" % (n,(doub+el.beta),el.doubFit[1]))
                    
                    
        if saveA == True:
            if case == 1:
                fig.savefig('Fig_4.1 Position and Intensity.png')
            if case == 2:
                fig.savefig('Fig_4.2 Power Law Fit.png')
            plt.close(fig)
    
    if testA == True:
        print("PartIV completed")
    return

###############################################################################################

def PartV(case,plotA,saveA,printA): # Absorption Edge
    if testA == True:
        print("PartV running")
    
    path1 = "C:\\Users\\Jacob\\Google Drive\\School\\Physics Degree\\Current Classes\\PHY 334 Advanced Lat 1\\X-Ray\\Excel Files\\Mo Thickness 1 to 6\\"
    path2 = "C:\\Users\\Jacob\\Google Drive\\School\\Physics Degree\\Current Classes\\PHY 334 Advanced Lat 1\\X-Ray\\Excel Files\\Carousel Thickness 4 CSV\\" 
    path3 = "C:\\Users\\Jacob/Google Drive\\School\\Physics Degree\\Current Classes\\PHY 334 Advanced Lat 1\\X-Ray\\Excel Files\\Carousel\\" 
    
    plt.close()    
    fig = plt.figure(case,figsize=fn.size)
    
    if case == 1: # 1) plot intensities for all thicknesses
        
        for t in fn.thickness:
            
            min_keV,max_keV = 5,25
            cps_min,cps_max = 0,1200
            LL,UL = fn.KeVToBin(fn.A,fn.B,min_keV),fn.KeVToBin(fn.A,fn.B,max_keV)
            xticks = [17.67,17.67+fn.Eq3(40,4)]
            
            data = fn.OpenCSV(path1+t.csv,LL+4,UL+4,int)
            X,Y = fn.XY(True,data,0,1)
            ax = fig.add_subplot(t.plot_params[0])
            ax.set_xlabel("keV",fontsize=fn.FS)
            ax.set_ylabel("CPS",fontsize=fn.FS)
            ax.set_xlim([min(X),max(X)])
            ax.set_ylim([cps_min,cps_max])
            
            ax.plot(X,Y,t.plot_params[2])
            ax.vlines(xticks, ymin=cps_min, ymax=cps_max, linewidth=.3)
            
            dL = fn.Leg(0,t.plot_params[2],t.name)
            ZrL = fn.Leg(0,'k','Zr doublet')
            fn.Legend(ax,[dL,ZrL],t.plot_params[1])
            
    if case == 2: # 2) find attenuation length by fitting Eq 4
        
        for t in fn.thickness:
            
            min_keV,max_keV = 5,25
            cps_min,cps_max = 0,1200
            LL,UL = fn.KeVToBin(fn.A,fn.B,min_keV),fn.KeVToBin(fn.A,fn.B,max_keV)
            
            data = fn.OpenCSV(path1+t.csv,LL+4,UL+4,int)
            X,Y = fn.XY(True,data,0,1)
            ax = fig.add_subplot(t.plot_params[0])
            ax.set_xlabel("keV",fontsize=fn.FS)
            ax.set_ylabel("CPS",fontsize=fn.FS)
            ax.set_xlim([min(X),max(X)])
            ax.set_ylim([cps_min,cps_max])
            
            rework = False
            if rework == True:
                I = Y[(np.abs(X-fn.Mo.alpha)).argmin()] # attenuated CPS
                lam_params = fn.LamFit(I,t.thick)
                t.lam = lam_params
                print(lam_params)
            
            ax.plot(X,Y,t.plot_params[2])
            ax.vlines([fn.Mo.alpha,fn.Mo.beta], ymin=cps_min, ymax=cps_max, linewidth=.3)
            ax.hlines(t.lam[2], xmin=min_keV, xmax=max_keV, color='k', linewidth=2)            
            
            dL = fn.Leg(0,t.plot_params[2],t.name)
            fitL = fn.Leg(0,'k','$\\lambda = %s$ inch' % t.lam[0])
            fn.Legend(ax,[dL,fitL],t.plot_params[1])

    if case == 3: # 3) find mu(E) from fitted attenuation length
        
        for el in fn.Wheel:
            
            LL,UL = 4,2051
            
            data1 = fn.OpenCSV(path3+el.csv,LL,UL,int)
            X1,Y1 = fn.XY(True,data1,0,1)
            
            data2 = fn.OpenCSV(path2+el.csv,LL,UL,int)
            X2,Y2 = fn.XY(True,data2,0,1)
            
            ax = fig.add_subplot(el.plot[0])
            ax.set_xlabel("keV",fontsize=fn.FS)
            ax.set_ylabel("CPS",fontsize=fn.FS)
            ax.plot(X1,Y1,'y')
            ax.plot(X2,Y2,el.plot[2])
            ax.set_xlim([min(X1),max(X1)])
            
            OL = fn.Leg(0,'y','Original')
            AL = fn.Leg(0,el.plot[2],'Attenuated')
            fn.Legend(ax,[OL,AL],2)
        
    if case == 4: #Attenuated spectrum of all samples from Zn w/ thickness = .04 inch
             
        X1,Y1 = fn.Edge(path3)
        X2,Y2 = fn.Edge(path2)
        
        X,Y = [],[]
        results = []
        for i in range(len(Y1)):
            y0,yA = Y1[i],Y2[i]
            R = ((y0-yA)/y0)
            PER = R * 100
            X.append(X1[i])
            Y.append(PER)
            results.append([X[i],Y1[i],Y2[i],PER])
        results = np.array(results)
            
        ax = fig.add_subplot(111)
        ax.set_xlabel("keV",fontsize=fn.FS)
        ax.set_ylabel("% CPS After Attenuation",fontsize=fn.FS)
        ax.set_xlim([min(X)*.9,max(X)*1.1])
        ax.set_ylim([0,max(Y)*1.1])
        ax.grid(True)
        
        EL = fn.Leg(0,'m','Absorption Edge ~ 17.4')
        fn.Legend(ax,[EL],1)
        
        printL = True
        if printL == True:
            print("[ keV , Original , Attenuated , percent remaining ] = %s" % results)
        
        ax.scatter(X,Y,color='g',marker='s',s=30)
        ax.vlines(17.443, ymin=0, ymax=62.1, color='m', linewidth=3)
        
            
    # 4) find the magnitude of the jump at the edge by assuming power law scaling for mu(E) above/blow edge, Eq 5
    if saveA == True:
        if case == 1:
            fig.savefig('Fig_5.1 Spectrum.png')
        if case == 2:
            fig.savefig('Fig_5.2 Attenuation Length.png')
        if case == 3:
            fig.savefig('Fig_5.3 Attenuation.png')
        if case == 4:
            fig.savefig('Fig_5.4 Absorption Edge.png')
        plt.close(fig)
    
    if printA == True:
        print("")
        print("PartV")
        
        if case == 1: 
            print("")
            print("Zn's K_beta is close to Mo's K_alpha. placing Zn in front of Mo attenuates the CPS")
            print("")
            
        if case == 2:
            print("Fitting Lambda")
            for t in fn.thickness:
                print("[ lambda , attenuation , model_CPS , error , error per] = %s" % t.lam)
            print("")
            
        if case == 4:
            print("Absorption Edge")
            print("[ keV , Original , Attenuated , percent Attenuated ] = %s" % results)
            print("")
        print("")
    
    if testA == True:
        print("PartV completed")
    return
       
