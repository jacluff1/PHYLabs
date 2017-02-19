import djak.gen as dg
import numpy as np 

# INSTRUMENT UNCERTAINTY

dA = .005 # V
dt = .5 # ms

def err4(name,t):
	"""calculates error for part4
	args
	----
	name: name that goes into returned djak.gen.var object
	t: djak.gen.var class object
	"""
	return dg.var(name,np.exp(t.val),np.exp(t.val)*t.err,'ms')

def err6(name,t):
	"""calculates error for part6
	args
	----
	name: name that goes into returned djak.gen.var object
	t: djak.gen.var class object
	"""
	return dg.var(name,(1/2)*t.val,(1/2)*t.err,'ms')

err7 = err6
 
err8a = err4

err8b = err6

def quad_err(name,alpha,beta,gamma,T_data,A_data):
	T = np.array(T_data)
	A = np.array(A_data)

	short = np.sqrt(beta.val**2 - (4*alpha.val*(gamma.val-np.average(A))))

	dt_q = dt/T

	one = dt_q

	two_1 = alpha.err/(2*alpha.val)
	two_2 = (1/alpha.val) * (beta.val + short)
	two_3 = (2*(gamma.val-A)) / short
	two = two_1 * (two_2 + two_3)

	three_1 = beta.err/(2*alpha.val)
	three_2 = 1 + (beta.val*beta.err)/short
	three = three_1 * three_2

	four = (gamma.val*gamma.err)/short

	five = (A*dA)/short

	dT1 = np.sqrt(one**2 + two**2 + three**2 + four**2 + five**2)
	tau1_err = np.exp(np.average(dT1)) * dt

	T1 = T + (1/(2*alpha.val)) * (beta.val - short)
	tau1_val = np.exp(np.average(T1))

	tau1 = dg.var(name,tau1_val,tau1_err,'ms')

	return tau1

def exp_err(tau_2,A_0,T_data,A_data,double='no'):
	T = np.array(T_data)
	A = np.array(A_data)
	
	one = dt/np.log(A_0.val/A)
	two = A_0.err/(A*np.log(A_0.val/A)**2)
	three = (A_0.val*dA)/(A**2 * np.log(A_0.val/A)**2)

	dt2 = np.average( np.sqrt(one**2 + two**2 + three**2) )
	t2 = np.average( T/np.log(A_0.val/A) )

	tau2_val = t2
	tau2_err = dt2
	if double != 'no':
		tau2_val *= 1/2
		tau2_err *= 1/2

	tau2 = dg.var(tau_2.name,tau2_val,tau2_err,'ms')

	return tau2





