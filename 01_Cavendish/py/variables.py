import Cav as cav
import numpy as np
from djak.gen import var
            
# angle calibration
S = var(np.array([43.6,28.4,18.3,11.9,7.8,5.1])/100, .5/100, 'meters')
L = var(762/100, .1, 'meters')
T = var([396., 385., 380., 380., 382., 382.5], 'std', 's') # Period
theta_s = var([1.88152248, 1.70676281, 1.73349876, 1.78764081], 'std', 'mRad') # displacment from equilibrium angle

# boom
l_b = var(14.9/100, .5/100, 'm') # length of boom
m_b = var(7.1/1000, .5/1000, 'kg') # mass of boom
w_b = var(1.5/100, .5/100, 'm') # width of boom

# cavendish box
w_t = var(.5/100, .5/100, 'm') # glass thickness
W = var(2.5/100, .5/100, 'm') # width of box

# m and m'
m_m = var(cav.m.m, .5/1000, 'kg') # mass of m
m_d = var(cav.m.d, .5/100, 'm') # diameter of m
m_r = var(m_d.val/2, m_d.err/2, 'm') # radius of m

# M and M'
M_m = var(cav.M.m, .5/1000, 'kg') # mass of M
M_d = var(cav.M.d, .5/100, 'm') # diameter of M
M_r = var(M_d.val/2, M_d.err/2, 'm') # radius of M

# fitting params
#theta_0, A, gamma, omega
par_m = var([0,10,1.5,(2*np.pi)/(.385)],None,None)
theta_e = var(0.0, 0.1, 'mRad') # equilibrium angle
A = var(10.5, 0.2, 'mRad') # initial amplitude of decaying wave

gamma = var(1.6, 0.1, '(10^3 s)^-1')
omega = var(16.3, 0.2, '(10^3 s)^-1') # angular frequency 

# misc
d = var((l_b.val-m_r.val)/2, cav.d_err(l_b,m_r), 'm') # distance from rotation axis to center of small sphere
I_b = (m_b.val/12)*(l_b.val**2 + w_b.val**2)
I_s = 2*(m_m.val*d.val**2 + (2/5)*m_m.val*m_r.val**2)
I = var(I_b + I_s, cav.I_err(m_m,d,m_r,m_b,l_b,w_b), 'kg m^2') # moment of inertia of boom, m and m'
beta = var(2*I.val*T.val, cav.beta_err(gamma,I), None) # decay parameter - 
R = var(W.val/2 + w_t.val + M_r.val, cav.R_err(W,w_t,M_r), 'm') # distance between center of Ms and m's
k = var(cav.K_s(I.val,T.val), cav.k_err(T,I), 'N m') # torsion constant
G_s = var(cav.G_s(theta_s.val), cav.G_s_err(k,theta_s,R,M_m,m_m,d), 'N m^2 kg^-2') # static gravitational constant 

# S,L,l_b,m_b,w_b,w_t,W,m_m,m_d,m_r,M_m,M_d,M_r,par_m,theta_e,A,gamma,omega,T,theta_s,d,I,beta,R,k,G_s