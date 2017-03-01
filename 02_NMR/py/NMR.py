import numpy as np 
import djak.plot as dp
import djak.gen as dg
import djak.fits as df 
import matplotlib.pyplot as plt
import err as err

# FITTING FUNCTIONS

def quadfit(par,t):
	return par[0]*(t-par[1])**2 + par[2]*(t-par[1]) + par[3]

def expfit(par,t):
	return par[0]*np.exp(-t/par[1]) + par[2]

def testfit(amount,index,par):
	par[index] = par[index] + amount
	return par

# EXPERIMENTS

def part1():
	return .138

def part3():
	a1,r1 = 3.5,100 # cm,ms
	a2,r2 = 1.9,10 # cm,ms
	x = np.array([r1,r2])
	x = np.log(x)
	y = np.array([a1,a2])

	X = np.linspace(1,100,1000)
	X = np.log(X)

	m = df.m_exp(x,y)
	b = df.b_exp(x,y)

	Y = m*X + b

	y_goal = a1/np.e 
	i_goal = dg.nearest(Y,y_goal)
	x_goal = X[i_goal]

	print("approx tau_1", np.exp(x_goal))

	fig = plt.figure(figsize=(15,7.5))
	plt.plot(x,y,'bs')
	plt.plot(X,Y,'-r')
	plt.plot(x_goal,y_goal,'g*')
	plt.show()
	return

def part4(plotA=True,printA=False,saveA=False,fs=20,size=(15,7.5)):
	#plt.xscale('log')

	tau = np.arange(7,26) # ms
	tau = np.log(tau)
	amp = np.array([3.8, 3.4, 2.9, 2.4, 2, 1.7, 1.2, .9, .5, .3, .4, .7, .9, 1.2, 1.5, 1.8, 2.0, 2.3, 2.6]) # V
	amp = amp[2:]
	tau = tau[2:]

	T = np.linspace(5,30,1000)
	T = np.log(T)
	par_m = [12.5,2.8,0,.5]
	A_m = quadfit(par_m,T)

	par_f = df.least_square(quadfit,par_m,tau,amp)
	# testfit(2,0,par_f) # alpha
	# testfit(.7,2,par_f) # beta
	# testfit(.2,3,par_f) # A_0
	# testfit(.05,1,par_f) # t_0
	A_f = quadfit(par_f,T)

	alpha = dg.var(par_f[0],2,'V/s^2')
	beta = dg.var(par_f[2],.7,'V/s')
	t_1 = dg.var(par_f[1],.05,'ms')
	A_0 = dg.var(par_f[3],.2,'V')

	# tau_1 = t_1
	tau_1 = err.err4(t_1)
	# tau1 = err.quad_err(alpha,beta,A_0,T,A_f)

	if plotA == True:
		fig = plt.figure(figsize=size)
		plt.title('$\\tau_1$ Two Pulse, Zero Crossing : Mineral Oil', fontsize=fs+2)
		plt.xlabel('ln(delay time/ms)', fontsize=fs)
		plt.ylabel('amplitude [V]', fontsize=fs)

		plt.plot(tau,amp, 's', color='b', label='data')
		# plt.plot(T,A_m,'-', color='m', lw=2, label='manual fit')
		plt.plot(T,A_f,'-', color='r', lw=2, label='least square fit')
		note='$A = \\alpha (t_q-t_1)^2 + \\beta (t_q-t_1) + \gamma$\n$\\alpha = %s \pm %s\ [V/s^2]$\n$\\beta = %s \pm %s\ [V/s]$\n$\\gamma = %s \pm %s\ [V]$\n$ln(t_1/ms) = %s \pm %s$\n$\\tau_1 = %s \pm %s\ [ms]$' % (alpha.pval,alpha.perr,beta.pval,beta.perr,A_0.pval,A_0.perr,t_1.pval,t_1.perr,tau_1.pval,tau_1.perr)
		plt.annotate(note, xy=(2.5,4), color='r', fontsize=fs)
	
		plt.xlim([min(T),max(T)])
		plt.legend(loc='best', numpoints=1)
		plt.tight_layout()
	
		if saveA == True:
			fig.savefig('../graphs/fig01.png')
			plt.close(fig)
		else:
			plt.show()

	if printA == True:
		print("tau_1 = %s +- %s" % (tau_1.pval,tau_1.perr))
		# print("tau1 = %s +- %s" % (tau1.pval,tau1.perr))

	return tau_1

def part5(plotA=True,printA=False,saveA=False,fs=20,size=(15,7.5)):

	tau = np.arange(1,17) # ms
	tau *= 2
	amp = np.array([7.5, 6.8, 6.0, 5.3, 4.9, 4.2, 4.0, 3.6, 3.4, 3.0, 2.8, 2.6, 2.4, 2.2, 2.0, 1.8]) # V

	T = np.linspace(0,30,1000)
	par_m = [8,10,.01]
	A_m = expfit(par_m,T)

	par_f = df.least_square(expfit,par_m,tau,amp)
	# testfit(.2,0,par_f) # A_0
	# testfit(.3,1,par_f) # tau_2
	# testfit(.2,2,par_f) # C
	A_f = expfit(par_f,T)

	A_0 = dg.var(par_f[0],.2,'V')
	tau_2 = dg.var(par_f[1],.3,'ms')
	# tau2 = err.exp_err(A_0,T[1:],A_f[1:],double='yes')
	C = dg.var(par_f[2],.2,'V')

	if plotA == True:
		fig = plt.figure(figsize=size)
		plt.title('$\\tau_2$ Two-Pulse Spin Echo : Mineral Oil', fontsize=fs+2)
		plt.xlabel('delay time [ms]', fontsize=fs)
		plt.ylabel('amplitude [V]', fontsize=fs)
		
		plt.plot(tau,amp, 's', color='b', label='data')
		# plt.plot(T,A_m,'-', color='m', lw=2, label='manual fit')
		plt.plot(T,A_f,'-', color='r', lw=2, label='least square fit')
		note='$A = A_0 e^{-t/\\tau_2} + C$\n$A_0 = %s \pm %s\ [V]$\n$\\tau_2 = %s \pm %s\ [ms]$\n$C = %s \pm %s\ [V]$' % (A_0.pval,A_0.perr,tau_2.pval,tau_2.perr,C.pval,C.perr)
		plt.annotate(note, xy=(15,4), color='r', fontsize=fs)
		
		plt.xlim([min(T),max(T)])
		plt.legend(loc='best', numpoints=1)
		plt.tight_layout()

		if saveA == True:
			fig.savefig('../graphs/fig02.png')
			plt.close(fig)
		else:
			plt.show()

	if printA == True:
		print("tau_2 = %s +- %s" % (tau_2.pval,tau_2.perr))
		# print("tau2 = %s +- %s" % (tau2.pval,tau2.perr))

	return tau_2

def part6(plotA=True,printA=False,saveA=False,fs=20,size=(15,7.5)):

	# tau = np.arange(20,120,20) # ms
	# amp = np.array([4.2, 2.10, 1.09, .62, .34]) # V
	tau = np.arange(2,18,2) # ms
	amp = np.array([9.52, 8.02, 6.5, 5.3, 3.96, 3.06, 2.36, 1.70])

	T = np.linspace(0,40,1000)
	par_m = [12,10,-1]
	A_m = expfit(par_m,T)

	par_f = df.least_square(expfit,par_m,tau,amp)
	# testfit(.5,0,par_f) # A_0
	# testfit(.3,1,par_f) # t_2
	# testfit(.2,2,par_f) # C
	A_f = expfit(par_f,T)

	A_0 = dg.var(par_f[0],.5,'V')
	t_2 = dg.var(par_f[1],.3,'ms')
	C = dg.var(par_f[2],.2,'V')

	tau_2 = t_2
	# tau_2 = err.err6(t_2)
	# tau2 = err.exp_err(A_0,T[1:],A_f[1:],double='yes')

	if plotA == True:
		fig = plt.figure(figsize=size)
		plt.title('$\\tau_2$ Carr-Purcell : Mineral Oil', fontsize=fs+2)
		plt.xlabel('delay time [ms]', fontsize=fs)
		plt.ylabel('amplitude [V]', fontsize=fs)

		plt.plot(tau,amp, 's', color='b', label='data')
		# plt.plot(T,A_m,'-', color='m', lw=2, label='manual fit')
		plt.plot(T,A_f,'-', color='r', lw=2, label='least square fit')
		note='$A = A_0 e^{-t/\\tau_2} + C$\n$A_0 = %s \pm %s\ [V]$\n$\\tau_2 = %s \pm %s\ [ms]$\n$C = %s \pm %s\ [V]$' % (A_0.pval,A_0.err,tau_2.pval,tau_2.perr,C.pval,C.perr)
		plt.annotate(note, xy=(20,3), color='r', fontsize=fs)

		plt.xlim([min(T),max(T)])
		plt.legend(loc='best', numpoints=1)
		plt.tight_layout()
	
		if saveA == True:
			fig.savefig('../graphs/fig03.png')
			plt.close(fig)
		else:
			plt.show()

	if printA == True:
		print("tau2 = %s +- %s" % (tau_2.pval,tau_2.perr))
		# print("tau2 = %s +- %s" % (tau2.pval,tau2.perr))

	return tau_2

def part7(plotA=True,printA=False,saveA=False,fs=20,size=(15,7.5)):

	# tau = np.arange(20,124,20) # ms
	# amp = np.array([4.20, 2.39, 1.37, 1.00, .65, .35]) # V
	tau = np.arange(2,32,2)
	amp = np.array([9.84, 8.86, 7.92, 7.06, 6.58, 6.06, 5.50, 5.12, 4.74, 4.50, 4.14, 3.96, 3.70, 3.50, 3.30])

	T = np.linspace(0,60,1000)
	par_m = [12,20,1]
	A_m = expfit(par_m,T)

	par_f = df.least_square(expfit,par_m,tau,amp)
	# testfit(.3,0,par_f) # A_0
	# testfit(.3,1,par_f) # t_2
	# testfit(.1,2,par_f) # C
	A_f = expfit(par_f,T)

	A_0 = dg.var(par_f[0],.3,'V')
	t_2 = dg.var(par_f[1],.3,'ms')
	C = dg.var(par_f[2],.1,'V')

	tau_2 = t_2
	# tau_2 = err.err7(t_2)
	# tau2 = err.exp_err(tau_2,A_0,T[1:],A_f[1:],double='yes')

	if plotA == True:
		fig = plt.figure(figsize=size)
		plt.title('$\\tau_2$ Meiboom-Gill : Mineral Oil', fontsize=fs+2)
		plt.xlabel('delay time [ms]', fontsize=fs)
		plt.ylabel('amplitude [V]', fontsize=fs)
	
		plt.plot(tau,amp, 's', color='b', label='data')
		# plt.plot(T,A_m,'-', color='m', lw=2, label='manual fit')
		plt.plot(T,A_f,'-', color='r', lw=2, label='least square fit')
		note='$A = A_0 e^{-t/\\tau_2} + C$\n$A_0 = %s \pm %s\ [V]$\n$\\tau_2 = %s \pm %s\ [ms]$\n$C = %s \pm %s\ [V]$' % (A_0.pval,A_0.perr,tau_2.pval,tau_2.perr,C.pval,C.perr)
		plt.annotate(note, xy=(30,4), color='r', fontsize=fs)

		plt.xlim([min(T),max(T)])
		plt.legend(loc='best', numpoints=1)
		plt.tight_layout()

		if saveA == True:
			fig.savefig('../graphs/fig04.png')
			plt.close(fig)
		else:
			plt.show()

	if printA == True:
		print("tau2 = %s +- %s" % (tau_2.pval,tau_2.perr))
		# print("tau2 = %s +- %s" % (tau2.pval,tau2.perr))

	return tau_2

def part8a(plotA=True,printA=False,saveA=False,fs=20,size=(15,7.5)):

	tau = np.arange(12,25) # ms
	tau = np.log(tau)
	amp = np.array([2.98, 2.56, 2.14, 1.78, 1.42, 1.09, .62, .59, .82, 1.07, 1.29, 1.59, 1.90]) # V

	T = np.linspace(10,25,1000)
	T = np.log(T)
	par_m = [10,np.log(19),0,.8]
	A_m = quadfit(par_m,T)

	par_f = df.least_square(quadfit,par_m,tau,amp)
	# testfit(1,0,par_f) # alpha
	# testfit(.7,2,par_f) # beta
	# testfit(.05,3,par_f) # A_0
	# testfit(.05,1,par_f) # t_0
	A_f = quadfit(par_f,T)

	alpha = dg.var(par_f[0],1,'V/s^2')
	beta = dg.var(par_f[2],.7,'V/s')
	t_1 = dg.var(par_f[1],.05,'ms')
	A_0 = dg.var(par_f[3],.05,'V')

	tau_1 = err.err8a(t_1)
	# tau1 = err.quad_err('part8a',alpha,beta,A_0,T,A_f)

	if plotA == True:
		fig = plt.figure(figsize=size)
		plt.title('$\\tau_1$ Two Pulse, Zero Crossing : Glycerin', fontsize=fs+2)
		plt.xlabel('ln(delay time/ms)', fontsize=fs)
		plt.ylabel('amplitude [V]', fontsize=fs)

		plt.plot(tau,amp, 's', color='b', label='data')
		# plt.plot(T,A_m,'-', color='m', lw=2, label='manual fit')
		plt.plot(T,A_f,'-', color='r', lw=2, label='least square fit')
		note='$A = \\alpha (t_q-t_1)^2 + \\beta (t_q-t_1) + \gamma$\n$\\alpha = %s \pm %s\ [V/s^2]$\n$\\beta = %s \pm %s\ [V/s]$\n$\\gamma = %s \pm %s\ [V]$\n$ln(t_1/ms) = %s \pm %s$\n$\\tau_1 = %s \pm %s\ [ms]$' % (alpha.pval,alpha.perr,beta.pval,beta.perr,A_0.pval,A_0.perr,t_1.pval,t_1.perr,tau_1.pval,tau_1.perr)
		plt.annotate(note, xy=(2.8,1.5), color='r', fontsize=fs)

		plt.xlim([min(T),max(T)]) 
		plt.legend(loc='best', numpoints=1)
		plt.tight_layout()

		if saveA == True:
			fig.savefig('../graphs/fig05.png')
			plt.close(fig)
		else:
			plt.show()

	if printA == True:
		print("tau_1_g = %s +- %s" % (tau_1.pval,tau_1.perr))
		# print("tau1_g = %s +- %s" % (tau1.pval,tau1.perr))

	return tau_1

def part8b(plotA=True,printA=False,saveA=False,fs=20,size=(15,7.5)):

	tau = np.arange(20,100,20) # ms
	amp = np.array([4, 1.71, .76, .43]) # V

	T = np.linspace(0,120,1000)
	par_m = [7,30,0]
	A_m = expfit(par_m,T)

	par_f = df.least_square(expfit,par_m,tau,amp)
	# testfit(.3,0,par_f) # A_0
	# testfit(.9,1,par_f) # t_2
	# testfit(.1,2,par_f) # C
	A_f = expfit(par_f,T)

	A_0 = dg.var(par_f[0],.3,'V')
	t_2 = dg.var(par_f[1],.9,'ms')
	C = dg.var(par_f[2],.1,'V')

	tau_2 = t_2
	# tau_2 = err.err8b(t_2)
	# tau2 = err.exp_err(A_0,T[1:],A_f[1:],double='yes')

	if plotA == True:
		fig = plt.figure(figsize=size)
		plt.title('$\\tau_2$ Meiboom-Gill : Glycerin', fontsize=fs+2)
		plt.xlabel('delay time [ms]', fontsize=fs)
		plt.ylabel('amplitude [V]', fontsize=fs)

		plt.plot(tau,amp, 's', color='b', label='data')
		# plt.plot(T,A_m,'-', color='m', lw=2, label='manual fit')
		plt.plot(T,A_f,'-', color='r', lw=2, label='least square fit')
		note='$A = A_0 e^{-t/\\tau_2} + C$\n$A_0 = %s \pm %s\ [V]$\n$\\tau_2 = %s \pm %s\ [ms]$\n$C = %s \pm %s\ [V]$' % (A_0.pval,A_0.perr,tau_2.pval,tau_2.perr,C.pval,C.perr)
		plt.annotate(note, xy=(60,1.5), color='r', fontsize=fs)

		plt.xlim([min(T),max(T)])
		plt.legend(loc='best', numpoints=1)
		plt.tight_layout()
	
		if saveA == True:
			fig.savefig('../graphs/fig06.png')
			plt.close(fig)
		else:
			plt.show()

	if printA == True:
		print("tau_2_g = %s +- %s" % (tau_2.pval,tau_2.perr))
		# print("tau2_g = %s +- %s" % (tau2.pval,tau2.perr))

	return tau_2

def RUN():
	part4(saveA=True)
	part5(saveA=True)
	part6(saveA=True)
	part7(saveA=True)
	part8a(saveA=True)
	part8b(saveA=True)
	return

