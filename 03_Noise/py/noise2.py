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

    R               =   np.array([ 1,           1e1,        1e2,        1e3,        1e4,        1e5,        1e6     ])
    G2              =   np.array([ 10*10*20,    10*10*20,   10*10*20,   10*10*20,   10*10*10,   1*10*40,    1*10*30 ])
    VOUT            =   np.array([ .9727, .9722, .9715, .9718, .9743, .9773, .9747, .9755, 1.0005, 1.0025, .9998, 1.0010,
                                1.2441, 1.2492, 1.2468, 1.2435, .9323, .9319, .9355, .9337, .8912, .8882, .8863, .8877,
                                .9585, .9550, .9595, .9650 ]).reshape([7,4])
    VOUT            =   np.array([ np.average(VOUT[i]) for i in range(7) ])

    Y_data          =   ( VOUT * (1e6)**2 ) / ( gain(G2)**2)                    # micro V^2
    X_fit           =   np.linspace(0, 2e4, 1000)                             # Ohm
    cutoff          =   2                                                       # num
    Y_fit, m, b     =   df.lin_fit( R[:-cutoff] , Y_data[:-cutoff] , X_fit )    # micro V^2
    R2              =   df.R2(R[:-cutoff],Y_data[:-cutoff])                     # num

    amp             =   np.sqrt(b.val)
    amp_err         =   b.err/(2*np.sqrt(b.val))
    ampN            =   dg.var( amp,amp_err )                                   # micro V

    fig = plt.figure(figsize=figsize)
    plt.title('Part III: Amplifier Noise', fontsize=fs+2)
    plt.xlabel('R [$\\Omega$]', fontsize=fs)
    plt.ylabel('<Vj$^2$> [$\\mu$ V$^2$]', fontsize=fs)
    plt.xlim(min(X_fit),max(X_fit))
    plt.plot(R[:-cutoff],Y_data[:-cutoff], 'bo', markersize=15, label='data')
    plt.plot(X_fit,Y_fit, color='r', lw=lw, label='fit')

    note = '$Vj^2 = mR + b$\n$m = %s \pm %s\ \\mu\ V^2/\\Omega$\n$b = %s \pm %s\ \\mu\ V^2$\n$R^2 = %s$' % (m.pval,m.perr,b.pval,b.perr,round(R2,6) )
    plt.annotate(note, xy=(1e3,35), color='r', fontsize=fs)

    note2 = 'Pre-Amp Noise = $\sqrt{b} \pm \\Delta b/(2 \sqrt{b})$\n  $= %s \pm %s\ \\mu V$' % (ampN.pval,ampN.perr)
    plt.annotate(note2, xy=(1e3,30), color='m', fontsize=fs+2)

    plt.legend(loc='best', fontsize=fs)
    if saveA == True:
        fig.savefig('../graphs/fig02.png')
        plt.close(fig)
    else:
        plt.show(fig)

    if printA == True:
        print('amplifier noise')
        print('R = %s Ohm' % R)
        print('Vj^2 = %s micro V^2' % Y_data )
        print("The best fit was found by not considering the last %s data points" % cutoff )
        print("n = %s +/- %s micro V^2/Ohm" % (m.pval,m.perr) )
        print("b = %s +/- %s micro V^2/Ohm")
        print('amplifer noise = %s +/- %s micro V^2' % (ampN.pval,ampN.perr) )
        print('')

    # experimental k

    DELTA_f             =   np.array([ 355,         1077,       3554,       10774,      35543,      107740 ])   # Hz
    G2                  =   np.array([ 10*10*60,    10*10*40,   10*10*20,   10*10*10,   1*10*60,    1*10*40 ])
    VOUT                =   np.array([ .8,          1.0,        .9,         .7,         .8,         .9 ])       # V
    R                   =   100e3                                                                               # Ohm
    Y_data              =   ( VOUT * (1e9)**2 ) / ( gain(G2)**2 * 4 * R * T_room )                              # nV^2/Ohm/K
    X_fit               =   np.linspace(0,120000,1000)                                                          # Hz
    Y_fit, m, b         =   df.lin_fit( DELTA_f , Y_data , X_fit )
    R2                  =   df.R2(DELTA_f,Y_data)

    kval                =   m.val/(1e9)**2                                                                      # J/K
    kval_err            =   m.err/(1e9)**2
    k                   =   dg.var( kval,kval_err )

    fig = plt.figure(figsize=figsize)
    plt.title('Part III: Bolzmann Constant', fontsize=fs+2)
    plt.xlabel('$\\Delta$ f [Hz]', fontsize=fs)
    plt.ylabel('<Vj$^2$>/R/T [nV$^2$/$\\Omega$/K]', fontsize=fs)
    plt.xlim(min(X_fit),max(X_fit))
    plt.plot(DELTA_f,Y_data, 'bo', markersize=15, label='data')
    plt.plot(X_fit,Y_fit, color='r', lw=lw, label='fit')

    note = '$V^2 = m f + b$\n$m = %s \pm %s\ nV^2/ \\Omega/ K / Hz$\n$b = %s \pm %s\ nV^2/ \\Omega/ K$\n$R^2 = %s$' % (m.pval,m.perr,b.pval,b.perr,round(R2,3) )
    plt.annotate(note, xy=(1.5e3,1.2), color='r', fontsize=fs)

    note2 = '$k = (m \pm \\Delta m)/(1_{+9})^2$\n  $= %s \pm %s\ J/K$' % (k.pval,k.perr)
    plt.annotate(note2, xy=(1.5e3,1), color='m', fontsize=fs+2)

    plt.legend(loc='best', fontsize=fs)
    if saveA == True:
        fig.savefig('../graphs/fig03.png')
        plt.close(fig)
    else:
        plt.show(fig)

    if printA == True:
        print("Bolzmann constant")
        print("m = %s +/- %s 1e18 J/K" % (m.pval,m.perr) )
        print("k = %s +/- %s J/K" % (k.pval,k.perr) )

    return

def part4(printA=True,saveA=True,figsize=(15,15),fs=20,lw=2):
    print("")
    print("Part IV:")

    R   =   np.array([ 10,          10e3 ])                                     # Ohm
    G2  =   np.array([ 10*10*20,    10*10*10 ])

    # JOHNSON NOISE AT T_room with A_ext and B_ext
    VOUT            =   np.array([ np.average([.9831, .9843, .9841, .9820]),    np.average([.7650, .7653, .7669, .7662]) ]) # V^2
    Y_data          =   ( VOUT * (1e6)**2 ) / (gain(G2)**2)                                                                 # (micro V)^2
    X_fit           =   np.linspace(0, 11e3, 1000)                                                                          # Ohm
    Y_fit, m, b     =   df.lin_fit( R , Y_data , X_fit )
    R2              =   df.R2(R,Y_data)

    mval            =   round(m.val,4)
    merr            =   5e-4

    bval            =   round(b.val,1)
    berr            =   .5

    noise           =   round(np.sqrt(b.val),1)                                                                              # micro V^2
    noise_err       =   round(berr/(2*np.sqrt(b.val)),1)
    # ampN            =   dg.var( noise,noise_err )

    fig = plt.figure(figsize=figsize)
    plt.title('Part IV: Amplifier Noise from A$_{ext}$ and B$_{ext}$', fontsize=fs+2)
    plt.xlabel('R [$\\Omega$]', fontsize=fs)
    plt.ylabel('<Vj$^2$> [$\\mu$V$^2$]', fontsize=fs)
    plt.xlim(min(X_fit),max(X_fit))
    plt.plot(R,Y_data, 'bo', markersize=15, label='data')
    plt.plot(X_fit,Y_fit, color='r', lw=lw, label='fit')

    note = '$Vj^2 = mR + b$\n$m = %s \pm %s\ \\mu V^2/\\Omega$\n$b = %s \pm %s\ \\mu V^2$\n$R^2 = %s$' % (mval,merr,bval,berr,round(R2,4) )
    plt.annotate(note, xy=(1e2,19), color='r', fontsize=fs)

    note2 = 'Pre-Amp Noise = $\sqrt{b} \pm \\Delta b/(2 \sqrt{b})$\n  $= %s \pm %s\ \\mu V$' % (noise,noise_err)
    plt.annotate(note2, xy=(1e2,17), color='m', fontsize=fs+2 )

    plt.legend(loc='best', fontsize=fs)
    if saveA == True:
        fig.savefig('../graphs/fig04.png')
        plt.close()
    else:
        plt.show()

    if printA == True:
        print('amplifier noise from Aext and Bext')
        print('it was found that Cext did not measure at an correct value and was discarded.')
        print('m = %s +/- %s micro V^2/Hz' % (mval,merr) )
        print('b = %s +/- %s micro V^2' % (bval,berr) )
        print('amplifer noise = %s +/- %s micro V' % (noise,noise_err) )
        print('')

    # experimental k at T = 77 K

    VOUT2           =   np.array([ np.average([.9825, .9832, .9804, .9827]),    np.average([.3909, .3905, .3905, .3911]) ]) # V^2
    Y_data2         =   ( VOUT2 * (1e9)**2 ) / (gain(G2)**2 * 4 * 77 * 110961 )                                             # nV^2 / K / Hz
    Y_fit2, m, b    =   df.lin_fit( R , Y_data2 , X_fit )
    R22             =   df.R2(R,Y_data2)

    bval            =   round(b.val,2)
    berr            =   .05

    mval            =   round(m.val,6)
    merr           =   5e-6

    kval            =   round(m.val/(1e9)**2,24)                                                                             # J/K
    kerr            =   round(merr/(1e9)**2,24)
    # k               =   dg.var( kval,kval_err )

    fig = plt.figure(figsize=figsize)
    plt.title('Part IV: Boltzmann Constant at 77 K', fontsize=fs+2)
    plt.xlabel('R [$\\Omega$]', fontsize=fs)
    plt.ylabel('[nVj$^2$/K/Hz]', fontsize=fs)
    plt.xlim(min(X_fit),max(X_fit))
    plt.plot(R,Y_data2, 'bo', markersize=15, label='data')
    plt.plot(X_fit,Y_fit2, color='r', lw=lw, label='fit')

    note = '$nVj^2/T/\\Delta f = mR + b$\n$m = %s \pm %s\ nV^2/K/Hz/\\Omega$\n$b = %s \pm %s\ nV^2/K/Hz$\n$R^2 = %s$' % (mval,merr,bval,berr,round(R22,4) )
    plt.annotate(note, xy=(1e2,3e-1), color='r', fontsize=fs)

    note2 = '$k = (m \pm \\Delta m)/(1_{+9})^2$\n  $= %s \pm %s\ J/K$' % (kval,kerr)
    plt.annotate(note2, xy=(1e2, 2.8e-1), color='m', fontsize=fs+2)

    plt.legend(loc='best', fontsize=fs)
    if saveA == True:
        fig.savefig('../graphs/fig05.png')
        plt.close()
    else:
        plt.show()

    if printA == True:
        print('Bolzmann Constant at 77 K')
        print('it was found that Cext did not measure at an correct value and was discarded.')
        print('experimental k = %s +/- %s J/m' % (kval,kerr) )
        print('')

def part5(printA=True,saveA=True,figsize=(20,15),fs=20,lw=3):

    f1      =   np.array([  10,     30,     100,    300,    1000,   3000  ])        # Hz
    f2      =   np.array([  .33,    1,      3.3,    10,     33,     100 ]) * 1000   # Hz
    Delta_f =   np.array([  355,    1077,   3554,   10774,  35343,  107740 ])       # Hz

    csvs = [ 'one' , 'two', 'three' , 'four', 'five' , 'six' ]
    for i in csvs:
        A = dg.CSV(i,skip_row=1)
        # A[:,1] = np.log10(A[:,1])

    def plot_axis(i):
        yx  =   -30
        yn  =   -90
        v   =   -40
        f   =   -35
        d   =   -70

        data = np.load('../npy/%s.npy' % csvs[i])
        X = np.log10(data[:,0])
        Y = data[:,1]

        ax = plt.subplot(3,2,i+1)
        ax.set_title('f$_1$ = %s Hz, f$_2$ = %s kHz, $\\Delta$ f = %s Hz' % (f1[i],f2[i]/1000,Delta_f[i]), fontsize=fs+2)
        ax.set_xlabel('log(Hz)', fontsize=fs)
        ax.set_ylabel('(dB)', fontsize=fs)
        ax.set_xlim([min(X),max(X)])
        ax.set_ylim([yn,yx])

        ax.plot(X,Y,'b',lw=lw)
        ax.vlines(np.log10(f1[i]), ymin=yn, ymax=v, color='r', lw=lw)
        ax.vlines(np.log10(f2[i]), ymin=yn, ymax=v, color='r', lw=lw)
        ax.hlines(d, xmin=np.log10(f1[i]), xmax=np.log10(f2[i]), color='m', lw=lw)
        ax.annotate('f$_1$', xy=(np.log10(f1[i]),f), color='r', fontsize=fs)
        ax.annotate('f$_2$', xy=(np.log10(f2[i]),f), color='r', fontsize=fs)
        ax.annotate('$\\Delta$ f = %s Hz' % int(f2[i]-f1[i]), xy=(np.log10(f1[i])+.1,d+2), color='m', fontsize=fs)

        return ax

    fig = plt.figure(figsize=figsize)
    for i in range(6):
        plot_axis(i)
    plt.tight_layout()

    if saveA == True:
        fig.savefig('../graphs/fig06.png')
        plt.close(fig)
    else:
        plt.show()

    if printA == True:
        per_err = ((f2-f1)-Delta_f)/Delta_f * 100
        print("")
        print("Part V:")
        print("f1 = %s Hz" % f1 )
        print("f2 = %s Khz" % (f2/1000) )
        print("delta f = %s Hz" % Delta_f )
        print("f2-f1 = %s Hz" % (f2-f1) )
        print("percent error = %s" % per_err )
        print("percent error = %s +/- %s" % (np.average(per_err),np.std(per_err)) )

    return
