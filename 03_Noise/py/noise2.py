import numpy as np
import matplotlib.pyplot as plt
import djak.fits as df
import djak.gen as dg
import pdb

from importlib import reload
reload(dg)
reload(df)

#=================================
# Constants
#---------------------------------

G1 = 600
G3 = 1/np.sqrt(10)

k_default = 1.38064852e-23
R_default = 10e3 # Ohm
T_room = 298 # Kelvin

#=================================
# Auxillary Functions
#---------------------------------

def gain(G2):
    return G1 * G2 * G3

def out2rms(G2,vout):
    G = gain(G2)
    return np.sqrt(vout)/G

def v_out(G2,delta_f,k,R,T=T_room):
    G = gain(G2)
    return G * np.sqrt( 4 * k * T * R * delta_f )

def v_rms(delta_f,k,R,T=T_room):
    return np.sqrt( 4 * k * T * R * delta_f )

def v2_conv(G2,vout): # Johnson Noise
    G = gain(G2)
    return vout/G**2

def v2_calc(delta_f,k,R,T=T_room): # Johnson Noise
    return 4 * R * k * T * delta_f

def err(v_conv,v_calc):
    abs_err = abs(v_calc - v_conv)
    per_err = 100 * (1 - v_calc/v_conv)
    v_av = (v_conv + v_calc)/2
    v = dg.var(v_av,abs_err)
    return v, per_err

#=================================
# Parts
#---------------------------------

def part2(printA=True):
    print("")
    print("Part II")

    vout_m      =   1.37                                    # V
    G2          =   400
    vrms_conv   =   out2rms(G2,vout_m)                      # V
    FWHM        =   .6e-6                                   # s
    delta_f     =   1/FWHM                                  # Hz
    vout_c      =   v_out(G2,delta_f,k_default,R_default)   # V
    vrms_c      =   v_rms(delta_f,k_default,R_default)      # V
    vout, operr =   err(vout_c,vout_m)
    vrms, iperr =   err(vrms_c,vrms_conv)

    if printA == True:
        print("x-y mode from multiplier to channel 2 looks like a parabola")
        print("after averaging over a long time, distribution looks similar to Gaussian.")
        print("FWHM = .6 microseconds, delta_f = 1/FWHM = %s Hz" % delta_f)
        print("v_out (measured) = %s V" % vout_m)
        print("v_out (calculated) = %s V" % vout_c)
        print("v_out = %s +/- %s (%s percent) V" % (vout.pval,vout.perr,operr) )
        print("v_rms (converted) = %s V" % vrms_conv)
        print("v_rms (calculated) = %s V" % vrms_c)
        print("v_rms = %s +/- %s (%s percent) V" % (vrms.pval,vrms.perr,iperr) )

    return

def part3(printA=True,saveA=True,figsize=(15,15), fs=20, lw=2):
    print("")
    print("Part III")

    # compare

    G2              =   1*10*100
    R               =   10e3                # Ohm
    delta_f         =   110961              # Hz
    vout            =   .94                 # V
    vrms_conv       =   out2rms(G2,vout)    # V
    vrms_calc       =   v_rms(delta_f,k_default,R_default)
    vrms, perr      =   err(vrms_conv,vrms_calc)

    if printA == True:
        print("compare")
        print("v_out = ", vout)
        print("v_rms_conv = ", vrms_conv)
        print("v_rms_calc = ", vrms_calc)
        print("v_rms = %s +/- %s V (%s percent)" % (vrms.pval,vrms.perr,perr) )
        print("")

    # amplifier noise

    # old values
    # R               =   np.array([ 1,           1e2,        1e3,        1e5,        1e6 ])      # Ohm
    # G2              =   np.array([ 10*10*20,    10*10*20,   10*10*20,   1*10*40,    1*10*30 ])
    # VOUT            =   np.array([ .98,         1.01,       1.26,       .89,        .96 ])      # V

    R               =   np.array([ 1,           1e1,        1e2,        1e3,        1e4,        1e5,        1e6     ])
    G2              =   np.array([ 10*10*20,    10*10*20,   10*10*20,   10*10*20,   10*10*10,   1*10*40,    1*10*30 ])
    VOUT            =   np.array([ .9727, .9722, .9715, .9718, .9743, .9773, .9747, .9755, 1.0005, 1.0025, .9998, 1.0010,
                                1.2441, 1.2492, 1.2468, 1.2435, .9323, .9319, .9355, .9337, .8912, .8882, .8863, .8877,
                                .9585, .9550, .9595, .9650 ]).reshape([7,4])
    VOUT            =   np.array([ np.average(VOUT[i]) for i in range(7) ])

    V2              =   v2_conv(G2,VOUT) * 1e12                                                 # micro V^2
    X_fit           =   np.linspace(0, 1.1e6, 1000)                                             # Ohm
    Y_fit, m, ampN  =   df.lin_fit( R , V2 , X_fit )
    R2              =   df.R2(R,V2)

    fig = plt.figure(figsize=figsize)
    plt.title('Part III: Amplifier Noise', fontsize=fs+2)
    plt.xlabel('R [$\\Omega$]', fontsize=fs)
    plt.ylabel('<V$_J^2$> [$\\mu$V$^2$]', fontsize=fs)
    plt.xlim(min(X_fit),max(X_fit))
    plt.plot(R,V2, 'bo', markersize=15, label='data')
    plt.plot(X_fit,Y_fit, color='r', lw=lw, label='fit')
    note = '$V_J^2 = mR + A$\n$m = %s \pm %s\ \\mu V^2/\\Omega$\n$A = %s \pm %s\ \\mu V^2$\n$R^2 = %s$' % (m.pval,m.perr,ampN.pval,ampN.perr,round(R2,4) )
    plt.annotate(note, xy=(1e4,250), color='r', fontsize=fs)
    plt.legend(loc='best', fontsize=fs)
    if saveA == True:
        fig.savefig('../graphs/fig02.png')
        plt.close(fig)
    else:
        plt.show(fig)

    if printA == True:
        print('amplifier noise')
        print('R = %s Ohm' % R)
        print('V^2 = %s micro V^2' % V2 )
        print('amplifer noise = %s +/- %s micro V' % (ampN.pval,ampN.perr) )
        print('')

    # experimental k

    DELTA_f             =   np.array([ 355,         1077,       3554,       10774,      35543,      107740 ])   # Hz
    G2                  =   np.array([ 10*10*60,    10*10*40,   10*10*20,   10*10*10,   1*10*60,    1*10*40 ])
    VOUT                =   np.array([ .8,          1.0,        .9,         .7,         .8,         .9 ])       # V
    R                   =   100e3                                                                               # Ohm
    V2                  =   (v2_conv(G2, VOUT)/(4*R*T) )                                                        # craziness
    X_fit               =   np.linspace(0,120000,1000)                                                          # Hz
    Y_fit, m, b         =   df.lin_fit( DELTA_f , V2 , X_fit )
    k, perr             =   err(k_default,m.val)
    R2                  =   df.R2(DELTA_f,V2)

    fig = plt.figure(figsize=figsize)
    plt.title('Part III: Bolzmann Constant', fontsize=fs+2)
    plt.xlabel('$\\Delta$ f [Hz]', fontsize=fs)
    plt.ylabel('<V$_J^2$>/R/T [V$^2$ $\\Omega^{-1}$ K$^{-1}$]', fontsize=fs)
    plt.xlim(min(X_fit),max(X_fit))
    plt.plot(DELTA_f,V2, 'bo', markersize=15, label='data')
    plt.plot(X_fit,Y_fit, color='r', lw=lw, label='fit')
    note = '$V^2 = k f + N$\n$k = %s \pm %s\ V^2\ \\Omega^{-1}\ K^{-1}\ Hz^{-1}$\n$N = %s \pm %s\ V^2\ \\Omega^{-1}\ K^{-1}$\n$R^2 = %s$' % (m.val,m.err,b.pval,b.perr,round(R2,4) )
    plt.annotate(note, xy=(1.5e3,1.2e-18), color='r', fontsize=fs)
    plt.legend(loc='best', fontsize=fs)
    if saveA == True:
        fig.savefig('../graphs/fig03.png')
        plt.close(fig)
    else:
        plt.show(fig)

    if printA == True:
        print("Bolzmann constant")
        print("k = %s +/- %s SI" % (m.pval,m.perr) )
        print("error in k = %s SI (%s percent)" % (k.perr,perr) )
        print(m.val,m.err)

    return

def part4(printA=True,saveA=True,figsize=(15,15),fs=20,lw=2):
    print("")
    print("Part IV:")

    R   =   np.array([ 10,          10e3 ])                                     # Ohm
    G2  =   np.array([ 10*10*20,    10*10*10 ])

    # JOHNSON NOISE AT T_room with A_ext and B_ext
    V_meter         =   np.array([ np.average([.9831, .9843, .9841, .9820]),    np.average([.7650, .7653, .7669, .7662]) ]) # V^2
    V_J_square      =   V_meter/gain(G2)**2 * (1e6)**2                                                                      # (micro V)^2
    X_fit           =   np.linspace(0, 11e3, 1000)                                                                          # Ohm
    Y_fit, m, ampN  =   df.lin_fit( R , V_meter , X_fit )
    print(m.val,m.err,ampN.val,ampN.err)

    R2              =   df.R2(R,V_meter)

    fig = plt.figure(figsize=figsize)
    plt.title('Part IV: Amplifier Noise from A$_{ext}$ and B$_{ext}$', fontsize=fs+2)
    plt.xlabel('R [$\\Omega$]', fontsize=fs)
    plt.ylabel('<Vj$^2$> [$\\mu$V$^2$]', fontsize=fs)
    plt.xlim(min(X_fit),max(X_fit))
    plt.plot(R,V_meter, 'bo', markersize=15, label='data')
    plt.plot(X_fit,Y_fit, color='r', lw=lw, label='fit')
    note = '$Vj^2 = mR + A$\n$m = %s \pm %s\ \\mu V^2/\\Omega$\n$A = %s \pm %s\ \\mu V^2$\n$R^2 = %s$' % (m.pval,m.perr,ampN.pval,ampN.perr,round(R2,4) )
    plt.annotate(note, xy=(1e3,.8), color='r', fontsize=fs)
    plt.legend(loc='best', fontsize=fs)
    # pdb.set_trace()
    if saveA == True:
        fig.savefig('../graphs/fig04.png')
        plt.close()
    else:
        plt.show()

    if printA == True:
        print('amplifier noise from Aext and Bext')
        print('it was found that Cext did not measure at an correct value and was discarded.')
        print('R = %s Ohm' % R)
        print('VJ^2 = %s micro V^2' % V_J_square )
        print('amplifer noise = %s +/- %s micro V' % (ampN.pval,ampN.perr) )
        print('')

    # experimental k at T = 77 K

    V_meter2            =   np.array([ np.average([.9825, .9832, .9804, .9827]),    np.average([.3909, .3905, .3905, .3911]) ]) # V^2
    Y_data2             =   V_meter2 / (gain(G2)**2 * 4 * 77 * 110961 )                                                                       # V^2 / K / Hz
    # V_J_square2         =   V_neter/gain(G2)**2 * (1e6)**2                                                                      # (micro V)^2
    Y_fit2, k, ampN2    =   df.lin_fit( R , Y_data2 , X_fit )
    R22                 =   df.R2(R,Y_data2)

    fig = plt.figure(figsize=figsize)
    plt.title('Part IV: Boltzmann Constant at 77 K', fontsize=fs+2)
    plt.xlabel('[$\\Omega$]', fontsize=fs)
    plt.ylabel('[Vj$^2$ / K / Hz]', fontsize=fs)
    plt.xlim(min(X_fit),max(X_fit))
    plt.plot(R,Y_data2, 'bo', markersize=15, label='data')
    plt.plot(X_fit,Y_fit2, color='r', lw=lw, label='fit')
    note = '$Vj^2/T/\\Delta_{f} = kR + N$\n$k = %s \pm %s\ J/K$\n$N = %s \pm %s\ V^2/K/Hz$\n$R^2 = %s$' % (k.pval,k.perr,ampN2.pval,ampN2.perr,round(R22,4) )
    plt.annotate(note, xy=(1e3,3e-19), color='r', fontsize=fs)
    plt.legend(loc='best', fontsize=fs)
    if saveA == True:
        fig.savefig('../graphs/fig05.png')
        plt.close()
    else:
        plt.show()

    if printA == True:
        print('Bolzmann Constant at 77 K')
        print('it was found that Cext did not measure at an correct value and was discarded.')
        print('R = %s Ohm' % R)
        print('Ydata2 = %s V^2/K/Hz' % Y_data2 )
        print('amplifer noise = %s V' % (ampN.val/(4*77*110961)) )
        print('experimental k = %s' % k.val )
        print('')

def part5(printA=True,saveA=True,figsize=(15,15),fs=20,lw=2):

    f1      =   np.array([  10,     30,     100,    300,    1000,   3000  ])    # Hz
    f2      =   np.array([  .33,    1,      3.3,    10,     33,     100 ])      # kHz
    Delta_f =   np.array([  355,    1077,   3554,   10774,  35343,  107740 ])   # Hz

    csvs = [ 'one' , 'two', 'three' , 'four', 'five' , 'six' ]
    npys = []
    for i in csvs:
        A = dg.CSV(i,skip_row=1,saveA=False)
        A[:,1] = np.log10(A[:,1])
        npys.append(A)

    def plot_axis(i):

        X = npy[i][:,0]
        Y = npy[i][:,1]

        ax = plt.subplot(3,2,i)
        ax.set_title('f$_1$ = %s Hz, f$_2$ = %s kHz, $\\Delta_f$ = %s Hz' % (f1[i],f2[i],Delta_f[i]), fontsize=fs+2)
        ax.set_xlabel('Hz', fontsize=fs)
        ax_set_ylabel('log(dB)', fontsize=fs)
        ax_set_xlim([min(X),max(X)])
        ax.plot(X,Y,'b',lw=lw)

        return ax

    fig = plt.figure(figsize=figsize)
    for i in range(6):
        plot_axis(i)
    plt.tight_layout()

    if saveA == True:
        fig.savefig('../graphs/06.png')
        plt.close(fig)
    else:
        plt.show()

    if printA == True:
        print("")
        print("Part V:")
        print("Oh My")
