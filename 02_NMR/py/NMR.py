import numpy as np 
import djak.plot as dp
import djak.gen as dg
import djak.fits as df 
import matplotlib.pyplot as plt 

def quadfit(par,t):
	return par[0]*(t-par[1])**2 + par[2]*(t-par[1]) + par[3]

def expfit(par,t):
	return par[0]*np.exp(-t/par[1])

def part4(plotA=True,printA=False,saveA=False,fs=20,size=(15,7.5)):

	fig = plt.figure(figsize=size)
	plt.title('$\\tau_1$ Two Pulses, Zero Crossing : Mineral Oil', fontsize=fs+2)
	plt.xlabel('delay time [ms]', fontsize=fs)
	plt.ylabel('amplitude [V]', fontsize=fs)
	#plt.xscale('log')

	tau = np.arange(7,26) # ms
	amp = np.array([3.8, 3.4, 2.9, 2.4, 2, 1.7, 1.2, .9, .5, .3, .4, .7, .9, 1.2, 1.5, 1.8, 2.0, 2.3, 2.6]) # V
	plt.plot(tau,amp, 's', color='b', label='data')

	T = np.linspace(5,30,1000)
	par_m = [.05,16,0,.5]
	A_m = quadfit(par_m,T)
	plt.plot(T,A_m,'-', color='m', lw=2, label='manual fit')

	#print(len(tau),len(amp))
	par_f = df.least_square(quadfit,par_m,tau,amp)
	A_f = quadfit(par_f,T)
	plt.plot(T,A_f,'-', color='r', lw=2, label='least square fit')

	alpha = par_f[0]
	tau_1 = par_f[1]
	A_0 = par_f[3]
	note='$A = \\alpha (t-\\tau_1)^2 + A_0$\n$\\alpha = %s\ [V]$\n$\\tau_1 = %s\ [ms]$\n$A_0 = %s\ [V]$' % (alpha,tau_1,A_0)
	plt.annotate(note, xy=(12,4), color='r', fontsize=fs)

	plt.xlim([min(T),max(T)])
	plt.legend(loc='best', numpoints=1)
	plt.tight_layout()
	
	if saveA == True:
		fig.savefig('../graphs/fig01.png')
		plt.close(fig)
	else:
		plt.show()


	return

def part5(plotA=True,printA=False,saveA=False,fs=20,size=(15,7.5)):

	fig = plt.figure(figsize=size)
	plt.title('$\\tau_2$ Two-Pulse Spin Echo : Mineral Oil', fontsize=fs+2)
	plt.xlabel('delay time [ms]', fontsize=fs)
	plt.ylabel('amplitude [V]', fontsize=fs)
	#plt.xscale('log')

	tau = np.arange(1,17) # ms
	amp = np.array([7.5, 6.8, 6.0, 5.3, 4.9, 4.2, 4.0, 3.6, 3.4, 3.0, 2.8, 2.6, 2.4, 2.2, 2.0, 1.8]) # V

	plt.plot(tau,amp, 's', color='b', label='data')

	T = np.linspace(0,30,1000)
	par_m = [8,10]
	A_m = expfit(par_m,T)
	plt.plot(T,A_m,'-', color='m', lw=2, label='manual fit')

	par_f = df.least_square(expfit,par_m,tau,amp)
	A_f = expfit(par_f,T)
	plt.plot(T,A_f,'-', color='r', lw=2, label='least square fit')

	A_0 = par_f[0]
	tau_2 = par_f[1]
	note='$A = A_0 e^{-t/\\tau_2}$\n$A_0 = %s\ [V]$\n$\\tau_2 = %s\ [ms]$' % (A_0,tau_2)
	plt.annotate(note, xy=(12,3), color='r', fontsize=fs)

	plt.xlim([min(T),max(T)])
	plt.legend(loc='best', numpoints=1)
	plt.tight_layout()
	
	if saveA == True:
		fig.savefig('../graphs/fig02.png')
		plt.close(fig)
	else:
		plt.show()

	return

def part6(plotA=True,printA=False,saveA=False,fs=20,size=(15,7.5)):

	fig = plt.figure(figsize=size)
	plt.title('$\\tau_2$ Carr-Purcell : Mineral Oil', fontsize=fs+2)
	plt.xlabel('delay time [ms]', fontsize=fs)
	plt.ylabel('amplitude [V]', fontsize=fs)
	#plt.xscale('log')

	tau = np.arange(20,120,20) # ms
	amp = np.array([4.2, 2.10, 1.09, .62, .34]) # V

	plt.plot(tau,amp, 's', color='b', label='data')

	T = np.linspace(0,120,1000)
	par_m = [8,30]
	A_m = expfit(par_m,T)
	plt.plot(T,A_m,'-', color='m', lw=2, label='manual fit')

	par_f = df.least_square(expfit,par_m,tau,amp)
	A_f = expfit(par_f,T)
	plt.plot(T,A_f,'-', color='r', lw=2, label='least square fit')

	A_0 = par_f[0]
	tau_2 = par_f[1]
	note='$A = A_0 e^{-t/\\tau_2}$\n$A_0 = %s\ [V]$\n$\\tau_2 = %s\ [ms]$' % (A_0,tau_2)
	plt.annotate(note, xy=(60,3), color='r', fontsize=fs)

	plt.xlim([min(T),max(T)])
	plt.legend(loc='best', numpoints=1)
	plt.tight_layout()
	
	if saveA == True:
		fig.savefig('../graphs/fig03.png')
		plt.close(fig)
	else:
		plt.show()

	return

def part7(plotA=True,printA=False,saveA=False,fs=20,size=(15,7.5)):

	fig = plt.figure(figsize=size)
	plt.title('$\\tau_2$ Meiboom-Gill : Mineral Oil', fontsize=fs+2)
	plt.xlabel('delay time [ms]', fontsize=fs)
	plt.ylabel('amplitude [V]', fontsize=fs)
	#plt.xscale('log')

	tau = np.arange(20,124,20) # ms
	amp = np.array([4.20, 2.39, 1.37, 1.00, .65, .35]) # V

	plt.plot(tau,amp, 's', color='b', label='data')

	T = np.linspace(0,120,1000)
	par_m = [7,40]
	A_m = expfit(par_m,T)
	plt.plot(T,A_m,'-', color='m', lw=2, label='manual fit')

	par_f = df.least_square(expfit,par_m,tau,amp)
	A_f = expfit(par_f,T)
	plt.plot(T,A_f,'-', color='r', lw=2, label='least square fit')

	A_0 = par_f[0]
	tau_2 = par_f[1]
	note='$A = A_0 e^{-t/\\tau_2}$\n$A_0 = %s\ [V]$\n$\\tau_2 = %s\ [ms]$' % (A_0,tau_2)
	plt.annotate(note, xy=(60,2), color='r', fontsize=fs)

	plt.xlim([min(T),max(T)])
	plt.legend(loc='best', numpoints=1)
	plt.tight_layout()
	
	if saveA == True:
		fig.savefig('../graphs/fig04.png')
		plt.close(fig)
	else:
		plt.show()

	return

def part8a(plotA=True,printA=False,saveA=False,fs=20,size=(15,7.5)):

	fig = plt.figure(figsize=size)
	plt.title('$\\tau_1$ Two Pulses, Zero Crossing : Glycerin', fontsize=fs+2)
	plt.xlabel('delay time [ms]', fontsize=fs)
	plt.ylabel('amplitude [V]', fontsize=fs)
	#plt.xscale('log')

	tau = np.arange(12,25) # ms
	amp = np.array([2.98, 2.56, 2.14, 1.78, 1.42, 1.09, .62, .59, .82, 1.07, 1.29, 1.59, 1.90]) # V
	plt.plot(tau,amp, 's', color='b', label='data')

	T = np.linspace(10,25,1000)
	par_m = [.05,19,0,.8]
	A_m = quadfit(par_m,T)
	plt.plot(T,A_m,'-', color='m', lw=2, label='manual fit')

	par_f = df.least_square(quadfit,par_m,tau,amp)
	A_f = quadfit(par_f,T)
	plt.plot(T,A_f,'-', color='r', lw=2, label='least square fit')

	alpha = par_f[0]
	tau_1 = par_f[1]
	A_0 = par_f[3]
	note='$A = \\alpha (t-\\tau_1)^2 + A_0$\n$\\alpha = %s\ [V]$\n$\\tau_1 = %s\ [ms]$\n$A_0 = %s\ [V]$' % (alpha,tau_1,A_0)
	plt.annotate(note, xy=(17,1.5), color='r', fontsize=fs)

	plt.xlim([min(T),max(T)])
	plt.legend(loc='best', numpoints=1)
	plt.tight_layout()
	
	if saveA == True:
		fig.savefig('../graphs/fig05.png')
		plt.close(fig)
	else:
		plt.show()

	return

def part8b(plotA=True,printA=False,saveA=False,fs=20,size=(15,7.5)):

	fig = plt.figure(figsize=size)
	plt.title('$\\tau_2$ Meiboom-Gill : Glycerin', fontsize=fs+2)
	plt.xlabel('delay time [ms]', fontsize=fs)
	plt.ylabel('amplitude [V]', fontsize=fs)
	#plt.xscale('log')

	tau = np.arange(20,100,20) # ms
	amp = np.array([4, 1.71, .76, .43]) # V

	plt.plot(tau,amp, 's', color='b', label='data')

	T = np.linspace(0,120,1000)
	par_m = [7,30]
	A_m = expfit(par_m,T)
	plt.plot(T,A_m,'-', color='m', lw=2, label='manual fit')

	par_f = df.least_square(expfit,par_m,tau,amp)
	A_f = expfit(par_f,T)
	plt.plot(T,A_f,'-', color='r', lw=2, label='least square fit')

	A_0 = par_f[0]
	tau_2 = par_f[1]
	note='$A = A_0 e^{-t/\\tau_2}$\n$A_0 = %s\ [V]$\n$\\tau_2 = %s\ [ms]$' % (A_0,tau_2)
	plt.annotate(note, xy=(60,1.5), color='r', fontsize=fs)

	plt.xlim([min(T),max(T)])
	plt.legend(loc='best', numpoints=1)
	plt.tight_layout()
	
	if saveA == True:
		fig.savefig('../graphs/fig06.png')
		plt.close(fig)
	else:
		plt.show()

	return