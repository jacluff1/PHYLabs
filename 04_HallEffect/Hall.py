import numpy as np
import pandas as pd
import djak.gen as dg
import djak.fits as df
import matplotlib.pyplot as plt

#===============================================================================
""" Global Variables and Parameters """
#-------------------------------------------------------------------------------

# plot parameters
pp      = {'figsize':(15,15),
           'fs':20,
           'lw':2,
           'ms':20,
           'data_style':'bo',
           'fit_style':'-r'}

# units
units   = {'time':'ms',
           'Temp':'C',
           'SampVolt':'V',
           'HallVolt':'mV',
           'SampCurr':'mA',
           'MagCurr':'A',
           'MagField':'kG'}

# constants and lab parameters
C       = {'N':         600,            # number of coils
           'mu0':       1/1000,         # permiability from Gs -> kGs
           'k':         8.6173303e-5,   # boltzmann constant [eV/K]
           'A':         np.pi * .04**2, # area of solenoid [m^2]
           'L':         60,             # inductance [H]
           'T':         .1,             # kGs -> T conversion
           'K':         273.15,         # K = C + 273.15
           'l_samp':    2,              # length of sample [cm]
           'a_samp':    .1              # cross section of sample [cm^2]
           }

#===============================================================================
""" Data """
#-------------------------------------------------------------------------------

# csv file: x = I_b, y = B
one_a   = {'g':.5, 'csv':'one_a'}
one_c   = {'g':1., 'txt':['data/one_c.txt','data/one_c_swap.txt']}
one_d   = {'g': np.array([ 0.3 , 0.5 , 1.0 , 1.5 , 2.0 ]),
           'B': np.array([ 1.7 , 1.1 , 0.6 , 0.4 , 0.3 ]),
           'csv':['one_d1', 'one_d2', 'one_d3', 'one_d4', 'one_d5'],
           'I_B': 1.00}

three_b = {'g':1.0, 'txt':'data/three_b.txt', 'I_p':29.8}
three_c = {'g':1.0, 'txt':'data/three_c.txt', 'I_p':29.8, 'I_B':0}
three_d = {'g':1.0, 'txt':'data/three_d.txt', 'I_p':30.1, 'I_B':2}

four_b = {'g':1.0, 'txt':'data/four_b.txt', 'I_p':30.1}
four_c = {'g':1.0, 'txt':'data/four_c.txt', 'I_p':30.1, 'I_B':0}
four_d = {'g':1.0, 'txt':'data/four_d.txt', 'I_p':30.1, 'I_B':1.99}

five_b = {'g':1.0, 'txt':'data/five_b.txt', 'I_p':9.1}
five_c = {'g':1.0, 'txt':'data/five_c.txt', 'I_p':9.1, 'I_B':0}
five_d = {'g':1.0, 'txt':'data/five_d.txt', 'I_p':9.1, 'I_B':2}


#===============================================================================
""" Lab Functions """
#-------------------------------------------------------------------------------

def equation1(I,gap):
    " find B"
    return C['N'] * C['mu0'] * I / gap

def equation5(rho,m,x,b):
    val     = np.exp(m.val * x + b.val)
    err     = val * np.sqrt( (x*m.err)**2 + (b.err)**2 )
    return dg.var(val,err)

#===============================================================================
""" Part 1.c)
Measure B vs current, for I B = 0 to Max to 0. You may want a short dwell-time,
like 100 ms. Use voltage control (current dial at max). This will prevent a nasty
spark and possible damage when you unplug the cables. Then measure opposite
polariy, by swapping the banana plugs at the supply (with V=0!), and run I B = 0 to
Max to 0. Note that R sense captures the polarity reversal. Swap the wires back
when finished. """
#-------------------------------------------------------------------------------

def fig01(saveA=True):

    """put units in for axies"""

    data_a      = pd.read_csv('data/one_a.csv',names=['I_B','B'])
    data_c      = pd.read_table(one_c['txt'][0])
    data_c_swap = pd.read_table(one_c['txt'][1])

    Bn          = data_c_swap['MagField'].values

    fig         = plt.figure(figsize=pp['figsize'])

    ax1         = plt.subplot(111)
    ax1.set_ylabel("B [%s]" % units['MagField'], fontsize=pp['fs'])
    ax1.set_xlabel("I$_B$ [%s]" % units['MagCurr'], fontsize=pp['fs'])

    ax1.set_title('Hysteresis Curve', fontsize=pp['fs']+2)
    ax1.plot(data_c['MagCurr'],data_c['MagField'],color='b', lw=pp['lw'], label='+ polarity')
    ax1.plot(data_c_swap['MagCurr'],data_c_swap['MagField'],color='g', lw=pp['lw'], label='- polarity')
    ax1.vlines(x=0,ymin=Bn[-1],ymax=Bn[0], color='r', lw=pp['lw']+2, label='remnant B')
    ax1.legend(loc='best', fontsize=pp['fs'])

    note = 'Remnant B = %.3f kGs' % (Bn[0] - Bn[-1])
    ax1.annotate(note, xy=(.1,0), color='r', fontsize=pp['fs'])

    plt.tight_layout()

    if saveA:   fig.savefig('png/fig01.png')
    else:       plt.show()

    dic         =   {'a':data_a,    'c':data_c,     'c_swap':data_c_swap,
                     'B0':Bn[0],    'Bf':Bn[-1],    'deltaB':Bn[0]-Bn[-1]}
    return      pd.Series(dic)

#===============================================================================
""" Part 1.d)
Measure B vs gap size "g" for a fixed current, say I B = 1 Amp. A single polarity
will suffice. Be sure that V = 0 when you loosen the pole pieces to change the gap
(maybe unplug a cable) Use spacer blocks to fix the gap, and tighten the pole
pieces to prevent accidental clamping shut (damage to sensor and fingers) """

""" Part 6.a)
Calculate the expected magnetic field, knowing N = 600 turns for each coil. (Hint:
using Ampere's law, and assume  r very large). Estimate the flux leakage (as %)
by comparing the calculated and measured magnetic field."""
#-------------------------------------------------------------------------------

def fig02(saveA=True):

    # data
    X           = one_d['g']
    Y1          = one_d['B']
    Y           = Y1*X**2

    # fit and theory
    X_fit       = np.linspace(.2,2.5,1000)
    Y_T         = equation1(1,X_fit) * X_fit**2

    # find flux leak
    Y_T1        = equation1(1,X) * X**2
    Y_L         = (Y_T1 - Y) / Y_T1
    yl_av       = np.average(Y_L)
    yl_std      = np.std(Y_L)
    yl          = dg.var(yl_av,yl_std)

    FIT         = df.lin_fit(X,Y,X_fit)
    Y_fit       = FIT['Y_fit']
    m           = FIT['m']
    b           = FIT['b']
    R2          = FIT['R2']

    fig         = plt.figure(figsize=pp['figsize'])
    plt.title('Magnetic Field Strength vs. Magnet Gap', fontsize=pp['fs']+2)
    plt.xlabel('Gap [cm]', fontsize=pp['fs'])
    plt.ylabel('B Gap$^2$ [%s cm$^2$]' % units['MagField'], fontsize=pp['fs'])
    plt.xlim([min(X_fit),max(X_fit)])

    plt.plot(X,Y,pp['data_style'], markersize=pp['ms'], label='data' )
    plt.plot(X_fit,Y_fit,pp['fit_style'], lw=pp['lw'], label='data fit' )
    plt.plot(X_fit,Y_T,'g', lw=pp['lw'], label='theory')

    plt.legend(loc='best', fontsize=pp['fs'], numpoints=1)

    note        = 'y = m x + b\n\
    m = %s $\pm$ %s [kG cm]\n\
    b = %s $\pm$ %s [kG cm$^2$]\n\
    R$^2$ = %s'\
    % (m.pval,m.perr,b.pval,b.perr,"{0:.4f}".format(R2) )
    plt.annotate(note, xy=(.3,1.1), color='r', fontsize=pp['fs'])

    note1       = '100 (B$_T$ - B$_L$) / B$_T$ = %s $\pm$ %s percent' % (yl.pval,yl.perr)
    plt.annotate(note1, xy=(.3,1), color='g', fontsize=pp['fs'])

    if saveA:
        fig.savefig('png/fig02.png')
        plt.close()
    else:
        plt.show()

    dic         =   {'X':X, 'Y':Y, 'fit':FIT, 'theory':Y_T, 'leak':yl.val, 'leak_error':yl.err}
    return pd.Series(dic)

#===============================================================================
""" Part 6.b)
Estimate  r (dimensionless ratio / 0 for the pole pieces, knowing the inductance
L = 60 Henry for the two coils in series."""
#-------------------------------------------------------------------------------

def fig07(saveA=True):

    data    = pd.read_table(three_d['txt'])
    B       = data['MagField'].values * C['T']      # T
    I       = data['MagCurr'].values                # A

    X       = B                                     # T
    Y       = I * C['L'] / C['A']                   # T

    X_fit       = np.linspace(min(X),max(X),1000)
    FIT         = df.lin_fit(X,Y,X_fit)
    Y_fit       = FIT['Y_fit']
    m           = FIT['m']
    b           = FIT['b']
    R2          = FIT['R2']

    mu_rval     = m.val
    mu_rerr     = m.err
    mu_r        = dg.var(mu_rval,mu_rerr)

    fig     = plt.figure(figsize=pp['figsize'])
    plt.title('Finding $\mu_r$', fontsize=pp['fs']+2)
    plt.xlabel('B [ T ]', fontsize=pp['fs'])
    plt.ylabel('I L / A [ T ]', fontsize=pp['fs'])
    plt.xlim([min(X),max(X)])

    plt.plot(X,Y, color='c', lw=pp['lw'], label='p-Ge')
    plt.plot(X_fit,Y_fit,pp['fit_style'], lw=pp['lw'], label='fit')
    plt.legend(loc='best', fontsize=pp['fs'])

    note        = 'y = m x + b\n\
    m = $\mu_r$ = %s $\pm$ %s [ H/m ]\n\
    b = %s $\pm$ %s [ T ]\n\
    R$^2$ = %s'\
    % (m.pval,m.perr,b.pval,b.perr,"{0:.4f}".format(R2) )
    plt.annotate(note, xy=(.097,24500), color='r', fontsize=pp['fs'])

    if saveA:
        fig.savefig('png/fig07.png')
        plt.close()
    else:
        plt.show()

    dic     =   {'mu_val':mu_r.val, 'mu_err':mu_r.err,  'mu_r':mu_r}
    return pd.Series(dic)

#===============================================================================
""" Part 6.c)
Estimate the magnetic remnance of the pole pieces."""
#-------------------------------------------------------------------------------

#===============================================================================
""" Part 6.d)
Plot Hall voltage vs field for all three samples on a single plot."""
#-------------------------------------------------------------------------------

def fig03(saveA=True):

    data1   = pd.read_table(three_b['txt'])
    X1      = data1['MagField'].values
    Y1      = data1['HallVolt'].values

    data2   = pd.read_table(four_b['txt'])
    X2      = data2['MagField'].values
    Y2      = data2['HallVolt'].values

    data3   = pd.read_table(five_b['txt'])
    X3      = data3['MagField'].values
    Y3      = data3['HallVolt'].values

    fig     = plt.figure(figsize=pp['figsize'])
    plt.title('Hall Voltage vs. Magnetic Field', fontsize=pp['fs']+2)
    plt.xlabel('B [%s]' % units['MagField'], fontsize=pp['fs'])
    plt.ylabel('V$_{Hall}$ [%s]' % units['HallVolt'], fontsize=pp['fs'])

    plt.plot(X1,Y1, color='c', lw=pp['lw'], label='p-Ge')
    plt.plot(X2,Y2, color='g', lw=pp['lw'], label='n-Ge')
    plt.plot(X3,Y3, color='b', lw=pp['lw'], label='i-Ge')
    plt.legend(loc='best', fontsize=pp['fs'])

    if saveA:
        fig.savefig('png/fig03.png')
        plt.close()
    else:
        plt.show()

    dic     =   {'d1':data1, 'd2':data2, 'd3':data3}
    return pd.Series(dic)

#===============================================================================
""" Part 6.e)
Plot sample voltage vs T for all three samples on a single plot."""
#-------------------------------------------------------------------------------

def fig04(saveA=True):

    data1   = pd.read_table(three_c['txt'])
    X1      = data1['Temp'].values
    Y1      = data1['SampVolt'].values

    data2   = pd.read_table(four_c['txt'])
    X2      = data2['Temp'].values
    Y2      = data2['SampVolt'].values

    data3   = pd.read_table(five_c['txt'])
    X3      = data3['Temp'].values
    Y3      = data3['SampVolt'].values

    fig     = plt.figure(figsize=pp['figsize'])
    plt.title('Sample Voltage vs. Temp', fontsize=pp['fs']+2)
    plt.xlabel('T [%s]' % units['Temp'], fontsize=pp['fs'])
    plt.ylabel('V$_{Sample}$ [%s]' % units['SampVolt'], fontsize=pp['fs'])

    plt.plot(X1,Y1, color='c', lw=pp['lw'], label='p-Ge')
    plt.plot(X2,Y2, color='g', lw=pp['lw'], label='n-Ge')
    plt.plot(X3,Y3, color='b', lw=pp['lw'], label='i-Ge')
    plt.legend(loc='best', fontsize=pp['fs'])

    if saveA:
        fig.savefig('png/fig04.png')
        plt.close()
    else:
        plt.show()

    dic     =   {'d1':data1, 'd2':data2, 'd3':data3}
    return pd.Series(dic)

#===============================================================================
""" Part 6.f)
Make a linearized plot for i-Ge, using the expected functional behavior, and from
this find the band gap (in eV) for Ge."""
#-------------------------------------------------------------------------------

def fig05(saveA=True):

    data        = pd.read_table(five_c['txt'])
    T           = data['Temp'].values               # K
    Vs          = data['SampVolt'].values           # V
    Is          = data['SampCurr'].values / 1000    # A
    l           = 20e-3                             # m
    A           = 10 / (1000)**2                    # m^2
    sigma       = l * Is / A / Vs                   # 1/(Ohm m)

    # linearize
    X           = 1/(T + 273.15)                    # 1/K
    Y           = np.log(Vs)

    X_fit       = np.linspace(min(X),max(X),1000)
    FIT         = df.lin_fit(X,Y,X_fit)
    Y_fit       = FIT['Y_fit']
    m           = FIT['m']
    b           = FIT['b']
    R2          = FIT['R2']

    E_gval      = 2 * C['k'] * m.val
    E_gerr      = 2 * C['k'] * m.err
    E_g         = dg.var(E_gval,E_gerr)

    fig = plt.figure(figsize=pp['figsize'])
    plt.title('i-Ge Band Gap', fontsize=pp['fs']+2)
    plt.xlabel("1/T [ 1/K ]", fontsize=pp['fs'])
    plt.ylabel("ln ( V$_s$ [ V ] )", fontsize=pp['fs'])
    plt.xlim(min(X_fit),max(X_fit))

    plt.plot(X,Y,'b', lw=pp['lw'], label='data')
    plt.plot(X_fit,Y_fit,'r', lw=pp['lw'], label='fit')

    note = '    Y=m x + b\n\
    ln( V$_s$ ) = ( -E$_g$/2k ) ( 1/T ) + ln( A )\n\
    m = %s $\pm$ %s\n\
    b = %s $\pm$ %s\n\
    R$^2$ = %.4f\n\
    E$_g$ = 2 k m $\pm$ 2 k $\Delta$ m\n\
         = %s $\pm$ %s eV'\
    % (m.pval,m.perr,b.pval,b.perr,R2,E_g.pval,E_g.perr)
    plt.annotate(note, xy=(.00255,.9), color='r', fontsize=pp['fs'])

    plt.legend(loc='best', numpoints=1)

    if saveA:
        fig.savefig('png/fig05.png')
        plt.close()
    else:
        plt.show()

    dic     =   {'m_val':m.val,     'm_err':m.err,      'm':m,
                 'E_g_val':E_g.val, 'E_gerr':E_g.err,   'E_g':E_g}
    return dic

#===============================================================================
""" Part 6.g)
Plot Hall voltage vs T for all three samples on a single plot. Why does the p-Ge
data change sign?"""
#-------------------------------------------------------------------------------

def fig06(saveA=True):

    data1   = pd.read_table(three_d['txt'])
    X1      = data1['Temp'].values
    Y1      = data1['HallVolt'].values

    data2   = pd.read_table(four_d['txt'])
    X2      = data2['Temp'].values
    Y2      = data2['HallVolt'].values
    # filter out blip
    # i_prob  = dg.maxima(X2,Y2,45,5)
    # print(i_prob,X2[i_prob],Y2[i_prob])
    X2      = np.delete(X2, [373,374,375])
    Y2      = np.delete(Y2, [373,374,375])
    data3   = pd.read_table(five_d['txt'])
    X3      = data3['Temp'].values
    Y3      = data3['HallVolt'].values

    fig     = plt.figure(figsize=pp['figsize'])
    plt.title('Hall Voltage vs. Temp', fontsize=pp['fs']+2)
    plt.xlabel('T [%s]' % units['Temp'] , fontsize=pp['fs'])
    plt.ylabel('V$_{Hall}$ [%s]' % units['HallVolt'], fontsize=pp['fs'])

    plt.plot(X1,Y1, color='c', lw=pp['lw'], label='p-Ge')
    plt.plot(X2,Y2, color='g', lw=pp['lw'], label='n-Ge')
    plt.plot(X3,Y3, color='b', lw=pp['lw'], label='i-Ge')
    plt.legend(loc='best', fontsize=pp['fs'])

    if saveA:
        fig.savefig('png/fig06.png')
        plt.close()
    else:
        plt.show()

    dic     =   {'d1':data1, 'd2':data2, 'd3':data3}
    return pd.Series(dic)

#===============================================================================
""" Part 6.h)
(opt) Find the mobilities  p and  n from the magneto-resistance data."""
#-------------------------------------------------------------------------------

#===============================================================================
""" Part 6.i)
Find the resistivity (ohm-cm) for p-Ge at T = 300K from the I-V data and the
sample geometry. From this, estimate the doping density (#/cm 3 )."""
#-------------------------------------------------------------------------------

def fig08(saveA=True):

    # raw data
    data    = pd.read_table(three_c['txt'])
    T       = data['Temp'].values + C['K']                  # temperature [K]
    Vs      = data['SampVolt'].values                       # sample voltage [V]
    Is      = data['SampCurr'].values / 1000                # sample current [A]

    # form data into axis
    X       = 1000/T                                        # [K]
    Y       = ( C['l_samp'] * Vs ) / ( C['a_samp'] * Is)    # [ohm-cm]
    Y       = np.log(Y)

    # start of semi-conductor phase (increase temperature -> lower resistance)
    semi0   = 2.95
    i_semi0 = dg.nearest(X,semi0)

    # fit semi-conductor phase to line
    X_fit   = np.linspace(2.7,3.4,1000)
    FIT     = df.lin_fit(X[ X>semi0 ],Y[ X>semi0 ],X_fit)
    Y_fit   = FIT['Y_fit']
    m       = FIT['m']
    b       = FIT['b']
    R2      = FIT['R2']

    fig     = plt.figure(figsize=pp['figsize'])
    plt.title('Resistivity of p-Ge', fontsize=pp['fs']+2)
    plt.xlabel('1000/T [ 1/K ]', fontsize=pp['fs'])
    plt.ylabel('ln( $\\rho$ [ $\Omega$ cm ] )', fontsize=pp['fs'])
    plt.xlim([min(X),max(X_fit)])
    plt.ylim([min(Y),max(Y_fit)])

    plt.plot(X[i_semi0:],Y[i_semi0:], color='b', lw=pp['lw'], label='p-Ge data')
    plt.plot(X_fit,Y_fit, pp['fit_style'], lw=pp['lw'], label='semi-phase fit')
    plt.vlines(x=semi0, ymin=6.5, ymax=8, color='r', linewidth=1)
    plt.legend(loc='best', fontsize=pp['fs'])

    note = 'semi-conductor phase fit cut off = %.0f K' % (1000/semi0)
    plt.annotate(note, xy=(semi0+.005,6.9), color='r', fontsize=pp['fs']-4, rotation=-90)

    note1 = 'y = m x + b\n\
    m = %s $\pm$ %s [ $\Omega$ cm K ]\n\
    b = %s $\pm$ %s [ $\Omega$ cm ]\n\
    R$^2$ = %.4f'\
    % (m.pval,m.perr,b.pval,b.perr,R2)
    plt.annotate(note1, xy=(3,6.99), color='r',fontsize=pp['fs']-4)

    T_interest  =   300 # K
    x_interest  =   1000/T_interest
    i_interest  =   dg.nearest(X_fit,x_interest)

    rho_val     =   np.exp(Y_fit[i_interest])
    rho         =   equation5(rho_val,m,x_interest,b)
    note2 = '$\\rho$ ( T = %.0f K ) = %s $\pm$ %s $\Omega$ cm' % (T_interest,rho.pval,rho.perr)
    plt.annotate(note2, xy=(3.17,6.85), color='g', fontsize=pp['fs']-4, rotation=-41)
    plt.plot(X_fit[i_interest],Y_fit[i_interest],'og', ms=pp['ms'])

    if saveA:
        fig.savefig('png/fig08.png')
        plt.close()
    else:
        plt.show()

    dic     =   {'mval':m.val,  'merr':m.err,   'm':m,
                 'bval':b.val,  'berr':b.err,   'b':b,
                 'R2':R2}
    return pd.Series(dic)
