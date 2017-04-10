import numpy as np
import pandas as pd
import djak.gen as dg
import djak.fits as df
import matplotlib.pyplot as plt

# Global Variables and Parameters===============================================

# plot parameters
pp      = {'figsize':(15,15),
           'fs':20,
           'lw':2,
           'ms':20,
           'data_style':'bo',
           'fit_style':'-r'}
Bunits  = 'T'

# Data==========================================================================

# csv file: x = I_b, y = B
one_a   = {'g':.5, 'csv':'one_a'}
one_c   = {'g':1., 'csv':['one_c','one_c_swap']}
one_d   = {'g': np.array([ 0.3 , 0.5 , 1.0 , 1.5 , 2.0 ]),
           'B': np.array([ 1.7 , 1.1 , 0.6 , 0.4 , 0.3 ]),
           'csv':['one_d1', 'one_d2', 'one_d3', 'one_d4', 'one_d5'],
           'I_B': 1.00}

three_b = {'g':1.0, 'txt':'data/three_b.txt', 'I_p':29.8}
three_c = {'g':1.0, 'txt':'data/three_c.txt', 'I_p':29.8, 'I_B':0}
three_d = {'g':1.0, 'txt':'data/three_d.txt', 'I_p':29.8, 'I_B':2}

four_b = {'g':1.0, 'txt':'data/four_b.txt', 'I_p':29.8}
four_c = {'g':1.0, 'txt':'data/four_c.txt', 'I_p':29.8, 'I_B':0}
four_d = {'g':1.0, 'txt':'data/four_d.txt', 'I_p':29.8, 'I_B':2}

five_b = {'g':1.0, 'txt':'data/five_b.txt', 'I_p':29.8}
five_c = {'g':1.0, 'txt':'data/five_c.txt', 'I_p':29.8, 'I_B':0}
five_d = {'g':1.0, 'txt':'data/five_d.txt', 'I_p':29.8, 'I_B':2}

# Ge      =   {}

# Aux Functions=================================================================


# Lab Functions=================================================================



# Part 1.c)=====================================================================
def fig01(saveA=True):

    """put units in for axies"""

    data_a      = pd.read_csv('data/one_a.csv',names=['I_B','B'])
    data_c      = pd.read_csv('data/one_c.csv',names=['I_B','B'])
    data_c_swap = pd.read_csv('data/one_c_swap.csv',names=['I_B','B'])

    fig         = plt.figure(figsize=pp['figsize'])

    ax1         = plt.subplot(121)
    ax1.set_title('Positive Polarity', fontsize=pp['fs']+2)
    ax1.plot(data_c['I_B'],data_c['B'],'bo', markersize=pp['ms'])

    ax2         = plt.subplot(122)
    ax2.set_title('Switched Polarity', fontsize=pp['fs']+2)
    ax2.plot(data_c_swap['I_B'],data_c_swap['B'],'go', markersize=pp['ms'])

    ax          = [ax1,ax2]
    for a in ax:
        a.set_xlabel("B [%s]" % Bunits, fontsize=pp['fs'])
        a.set_ylabel("I$_B$ [A]", fontsize=pp['fs'])

    plt.tight_layout()

    if saveA:   fig.savefig('png/fig01.png')
    else:       plt.show()

# Part 1.d)=====================================================================
def fig02(saveA=True):

    X           = one_d['g']
    Y           = one_d['B']

    X_fit       = np.linspace(.2,2.5,1000)

    par1        = [2,.6,.2]
    def fitfunc1(p,x):
        return p[0] * np.exp(-x/p[1]) + p[2]
    par_fit1    = df.least_square(fitfunc1,par1,X,Y)
    Y_fit1      = fitfunc1(par_fit1,X_fit)

    par2        = [1,1,1]
    def fitfunc2(p,x):
        return p[0]*x**p[1] + p[2]
    par_fit2    = df.least_square(fitfunc2,par2,X,Y)
    Y_fit2      = fitfunc2(par_fit2,X_fit)

    fig         = plt.figure(figsize=pp['figsize'])
    plt.title('I$_B$ vs. Gap', fontsize=pp['fs']+2)
    plt.xlabel('Gap [cm]', fontsize=pp['fs'])
    plt.ylabel('I$_B$ [%s]' % Bunits, fontsize=pp['fs'])
    plt.xlim([min(X_fit),max(X_fit)])

    plt.plot(X,Y,pp['data_style'], markersize=pp['ms'], label='data' )
    plt.plot(X_fit,Y_fit1,pp['fit_style'], lw=pp['lw'], label='exponential fit' )
    plt.plot(X_fit,Y_fit2,'m', lw=pp['lw'], label='power fit' )

    plt.legend(loc='best', fontsize=pp['fs'], numpoints=1)

    note        = '$I_B = c_1 x^{c_2} + c_3$\n$c_1 = %s \pm %s$\n$c_2 = %s \pm %s$\n$c_3 = %s \pm %s$' % (par_fit2[0],0,par_fit2[1],0,par_fit2[2],0)
    plt.annotate(note, xy=(.5,1.5), color='m', fontsize=pp['fs'])

    if saveA:
        fig.savefig('png/fig02.png')
        plt.close()
    else:
        plt.show()


# Part 6.a)=====================================================================
"""Calculate the expected magnetic field, knowing N = 600 turns for each coil. (Hint:
using Ampere's law, and assume  r very large). Estimate the flux leakage (as %)
by comparing the calculated and measured magnetic field."""

# Part 6.b)=====================================================================
"""Estimate  r (dimensionless ratio / 0 for the pole pieces, knowing the inductance
L = 60 Henry for the two coils in series."""

# Part 6.c)=====================================================================
"""Estimate the magnetic remnance of the pole pieces."""

# Part 6.d)=====================================================================
"""Plot Hall voltage vs field for all three samples on a single plot."""

def fig03(saveA=True):

    data1   = pd.read_table(three_b['txt'])
    X1      = data1['MagField'].values
    Y1      = data1['HallVolt'].values

    # data2   = pd.read_table(four_b['txt'])
    # X2      = data2['MagField'].values
    # Y2      = data2['HallVolt'].values
    #
    # data3   = pd.read_table(five_b['txt'])
    # X3      = data3['MagField'].values
    # Y3      = data3['HallVolt'].values

    fig     = plt.figure(figsize=pp['figsize'])
    plt.title('Hall Voltage vs. Magnetic Field', fontsize=pp['fs']+2)
    plt.xlabel('B [%s]' % Bunits, fontsize=pp['fs'])
    plt.ylabel('V$_{Hall}$ [V]', fontsize=pp['fs'])

    plt.plot(X1,Y1,'co', markersize=pp['ms'], label='p-Ge')
    # plt.plot(X2,Y2,'go', markersize=pp['ms'], label='n-Ge')
    # plt.plot(X3,Y3, 'bo', markersize=pp['ms'], label='i-Ge')

    if saveA:
        fig.savefig('png/fig03.png')
        plt.close()
    else:
        plt.show()

    return

# Part 6.e)=====================================================================
"""Plot sample voltage vs T for all three samples on a single plot."""

def fig04(saveA=True):

    data1   = pd.read_table(three_c['txt'])
    X1      = data1['Temp'].values
    Y1      = data1['SampVolt'].values

    # data2   = pd.read_table(four_c['txt'])
    # X2      = data2['Temp'].values
    # Y2      = data2['SampVolt'].values
    #
    # data3   = pd.read_table(five_c['txt'])
    # X3      = data3['Temp'].values
    # Y3      = data3['SampVolt'].values

    fig     = plt.figure(figsize=pp['figsize'])
    plt.title('p-Ge: Sample Voltage vs. Temp', fontsize=pp['fs']+2)
    plt.xlabel('T [K]', fontsize=pp['fs'])
    plt.ylabel('V$_{Sample}$ [V]', fontsize=pp['fs'])

    plt.plot(X1,Y1,'co', markersize=pp['ms'], label='p-Ge')
    # plt.plot(X2,Y2,'go', markersize=pp['ms'], label='n-Ge')
    # plt.plot(X3,Y3,'bo', markersize=pp['ms'], label='i-Ge')

    if saveA:
        fig.savefig('png/fig04.png')
        plt.close()
    else:
        plt.show()

# Part 6.f)=====================================================================
"""Make a linearized plot for i-Ge, using the expected functional behavior, and from
this find the band gap (in eV) for Ge."""

def fig05(saveA=True):

    data        = NotImplemented
    X           = data[NotImplemented].values
    Y           = data[NotImplemented].values

    X_fit       = np.linspace(min(X),max(X),1000)
    FIT         = df.lin_fit(X,Y,X_fit)
    Y_fit       = FIT['Y_fit']
    m           = FIT['m']
    b           = FIT['b']
    R2          = FIT['R2']

    fig = plt.figure(figsize=pp['figsize'])
    plt.title('i-Ge Band Gap', fontsize=pp['fs']+2)
    plt.xlabel(NotImplemented, fontsize=pp['fs'])
    plt.ylabel(NotImplemented, fontsize=pp['fs'])

    plt.plot(X,Y,'bo', markersize=pp['ms'], label='data')
    plt.plot(X_fit,Y_fit,'r', lw=pp['lw'], label='fit')

    note='$Y=m x + b$\n$m = %s \pm %s\ units$\n$b = %s \pm %s\ units$\n$R^2 = %s$' % (m.pval,m.perr,b.pval,b.perr,R2)
    plt.annotate(note, xy=NotImplemented, color='r', fontsize=pp['fs'])

    plt.legend(loc='best', numpoints=1)

# Part 6.g)=====================================================================
"""Plot Hall voltage vs T for all three samples on a single plot. Why does the p-Ge
data change sign?"""

def fig06(saveA=True):

    # data1   = pd.read_table(three_d['txt'])
    # X1      = data1['Temp'].values
    # Y1      = data1['HallVolt'].values

    # data2   = pd.read_table(four_d['txt'])
    # X2      = data2['Temp'].values
    # Y2      = data2['HallVolt'].values
    #
    # data3   = pd.read_table(five_d['txt'])
    # X3      = data3['Temp'].values
    # Y3      = data3['HallVolt'].values

    fig     = plt.figure(figsize=pp['figsize'])
    plt.title('Hall Voltage vs. Temp', fontsize=pp['fs']+2)
    plt.xlabel('T [K]', fontsize=pp['fs'])
    plt.ylabel('V$_{Hall}$ [V]', fontsize=pp['fs'])

    # plt.plot(X1,Y1,'co', markersize=pp['ms'], label='p-Ge')
    # plt.plot(X2,Y2,'go', markersize=pp['ms'], label='n-Ge')
    # plt.plot(X3,Y3,'bo', markersize=pp['ms'], label='i-Ge')

    if saveA:
        fig.savefig('png/fig06.png')
        plt.close()
    else:
        plt.show()

# Part 6.h)=====================================================================
"""(opt) Find the mobilities  p and  n from the magneto-resistance data."""

# Part 6.i)=====================================================================
"""Find the resistivity (ohm-cm) for p-Ge at T = 300K from the I-V data and the
sample geometry. From this, estimate the doping density (#/cm 3 )."""
