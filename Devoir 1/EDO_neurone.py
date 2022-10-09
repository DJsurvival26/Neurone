from scipy.integrate import solve_ivp as ode45 
from matplotlib import pyplot as plt
import numpy as np
import time
import matplotlib.pyplot as plt
import const as c

global I_app
I_app=100
global V_0
V_0=-40
global N_0
N_0=0.1
def M_inf(V):
    a=(1+np.tanh((V-c.V1)/c.V2))*1/2
    return a

def N_inf(V):
    a=(1+np.tanh((V-c.V3)/c.V4))*1/2
    return a

def Tau_N(V):
    a=1/(c.phi*(np.cosh((V-c.V3)/(2*c.V4))))
    return a

def dx(t,x):
    V=x[0]
    N=x[1]
    dx=np.zeros(2)
    dx[0]=((-c.gl*(V-c.El)-c.g_Ca*M_inf(V)*(V-c.E_Ca)-c.gk*N*(V-c.Ek)+I_app))*1/c.C
    dx[1]=(N_inf(V)-N)/Tau_N(V)
    print(dx)
    return dx

time=np.linspace(0,500,20000)
a=ode45(dx,[0,500],[V_0,N_0],t_eval=time)
V=a.y[0,:]
N=a.y[1,:]

plt.plot(time,V)
plt.plot(time,N)
plt.show



            
    