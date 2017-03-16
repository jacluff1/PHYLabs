import numpy as np
import matplotlib.pyplot as plt
import djak.fits as df
import djak.gen as dg

from importlib import reload
reload(dg)
reload(df)



#=================================
# Constants

G1 = 600
G3 = 1/np.sqrt(10)

k_default = 1.38064852e-23
R_default = 10e3 # Ohm
T = 295 # Kelvin

#---------------------------------



#=================================
# Auxillary Functions

def gain(G2):
    return G1 * G2

def V2(V_out):
    return np.sqrt(V_out*10)

def V1(G2,delta_f,k,R):
    G = gain(G2)
    return G*np.sqrt(4*k*T*R*delta_f)

def err(V1,V2):
    abs_err = abs(V1-V2)
    per_err = 100 * (1 - V2/V1)
    V_av = (V1 + V2)/2
    V = dg.var(V_av,abs_err)
    return V, per_err
#---------------------------------



#=================================
# Parts

def part2():
    print("")
    print("Part II")

    # section 2
    print("")
    v_out = 1.37
    G2 = 400
    result2 = V2(1.37)
    print("II b) v_rms = V2(v_out = %s ) = %s V" % (v_out,result2) )

    # section 3
    print("")
    print("II c) x-y mode from multiplier to channel 2 looks like a parabola")

    # section 4
    print("")
    FWHM = .6e-6 # seconds
    result4 = 1/FWHM # Hz
    print("II d) after averaging over a long time, distribution looks similar to Gaussian.")
    print("FWHM = .6 microseconds, delta_f = 1/FWHM = %s Hz" % result4)

    # section 5
    print("")
    G2 = 400
    result5 = V1(G2,result4,k_default,R_default) # V
    print("II e) v_rms = V1( G2 = %s , delta_f = %s Hz , k,R = default SI ) = %s V" % (G2,result4,result5) )

    print("")
    V, per_err = err(result5,result2)
    print("percent error = ", per_err)
    print("V_rms = %s +/- %s V" % (V.pval,V.perr) )

    return V, per_err

def part3(figsize=(15,15), fs=20, lw=2):
    print("")
    print("Part III")

    # section 4 a
    print("")
    f1 = 100 # Hz
    f2 = 100e3 # Hz
    delta_f = 110961 # Hz
    G2 = 1*10*100
    R = 10e3 # Ohm
    V_1 = V1(G2,delta_f,k_default,R_default)
    V_2 = V2(delta_f)
    V, per_err = err(V_1,V_2)
    print("V1 = ", V_1)
    print("V2 = ", V_2)
    print("per_err = ", per_err)
    print("V = %s +/- %s V" % (V.pval,V.perr) )




    # # section 4 b
    # print("")
    # R = np.array([1 , 1e2 , 1e3 , 1e5 , 1e6 ]) # Ohm
    # G2 = np.array([10*10*20 , 10*10*20 , 10*10*20 , 1*10*40 , 1*10*30 ])
    # V_2 = np.array([ V2(g2,delta_f,R=r) for g2,r in zip(G2,R) ])
    # print("III 4 b) V2^2 = %s" % (V_2)**2)
    # print("see fig01")
    #
    # fig = plt.figure(figsize=figsize)
    # plt.xlabel("Resistance [$\\Omega$]", fontsize=fs)
    # plt.ylabel("Voltage$^2$ [V$^2$]", fontsize=fs)
    # plt.plot(R,V_2**2, 'bo', markersize=15, label='data')
    #
    # X_fit = np.linspace(0,1100000,1000)
    # Y_fit , m , b = df.lin_fit( R, (V_2**2), X_fit )
    # plt.plot(X_fit,Y_fit, color='r', lw=lw, label='linear fit')
    # plt.xlim(min(X_fit),max(X_fit))
    #
    # note='$V^2 = m R + N$\n$m = %s \pm %s\ units$\n$N = %s \pm %s units$' % (m.pval,m.perr,b.pval,b.perr)
    # plt.annotate(note, xy=(600000,2), color='r', fontsize=fs)
    #
    # plt.legend(loc='best', fontsize=fs, numpoints=1)
    # plt.show()
