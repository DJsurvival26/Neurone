
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

V1 = -1.2
V2 = 18
V3 = 2
V4 = 30

phi = 0.04

C = 20

gk = 8
gl = 2
gca = 4

Eca = 120
Ek = -80
El = -60

def Minf (V) :
    return 0.5*(1+np.tanh((V-V1)/(V2)))

def Ninf (V) :
    return 0.5*(1+np.tanh((V-V3)/V4))

def tauN (V) :
    return (1/phi)*(1/(np.cosh((V-V3)/(2*V4))))

def syst (t, y, Iapp) : 
    dy = np.zeros(2)
    
    dy[0] = (1/C)*(-gl*(y[0]-El)-gca*Minf(y[0])*(y[0]-Eca)-gk*y[1]*(y[0]-Ek) + Iapp)
    dy[1] = (1/tauN(y[0]))*(Ninf(y[0])-y[1])
    
    return dy

def graphe_sol (Iapp, y0) :
    t_span = np.array([0, 500])
    
    sol = solve_ivp(lambda t, y: syst(t, y, Iapp), t_span, y0, max_step = 1)
    
    plt.figure(figsize = (15, 10))
    plt.xlim(0, 500)
    plt.ylim(-65, 40)
    plt.title('Solution de V au cours du temps avec Iapp = ' + str(Iapp))
    plt.xlabel('Temps (ms)')
    plt.ylabel('tension (mV)')
    
    plt.plot(sol.t, sol.y[0])
    print(sol.y[0][-1])
    plt.show()


def Null1 (V, Iapp) :
    return Ninf(V)
def Null2 (V, Iapp) :
    return (1/((V-Ek))*gk)*(-gl*(V-El)-gca*Minf(V)*(V-Eca) + Iapp)

def plan_de_phase (Iapp, y0) :
    V = np.linspace(-100, 100, 500)
    
    N1 = Null1(V, Iapp)
    N2 = Null2(V, Iapp)
    
    t_span = np.array([0, 500])
    
    sol = solve_ivp(lambda t, y: syst(t, y, Iapp), t_span, y0, max_step = 1)
    
    plt.figure(figsize = (15, 10))
    plt.xlim(-100, 100)
    plt.ylim(0, 1)
    
    plt.title('Plan de phase avec Iapp = ' + str(Iapp))
    plt.xlabel('V (mV)')
    plt.ylabel('N (ad)')
    
    plt.plot(V, N1,"red")
    plt.plot(V, N2,"yellow")
    plt.plot(sol.y[0], sol.y[1],"purple")
    
    plt.show
plt.xlim(-100, 100)
plt.ylim(0, 1)
plan_de_phase(100, [-40,0.1])
    
print(Null2(1.5,100))