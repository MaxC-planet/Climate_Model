#variables in governing equation 

import numpy as np
import matplotlib.pyplot as plt
import scipy

time=[i for i in range(365)]
T=[220+i for i in range(100)]

#lam: latitude

def S_11(lam, t, a, e, phi_ini, m_star, L): #stellar flux for tidally locked
    period=np.sqrt((4*(np.pi**2)*((a*(1.49598*(10**11)))**3))/((6.67259*(10**-11))*m_star))
    E_ini=2*np.arctan(np.sqrt(((1-e)/(1+e)))*np.tan(phi_ini/2))
    M_ini=E_ini-e*np.sin(E_ini)
    d_M=2*np.pi*t/period
    M=M_ini+d_M
    E=M
    for i in range(8):
        E=M+e*np.sin(E)
    r=a*(1-e*np.cos(E))
    q_0=L/(4*np.pi*((1.49598*(10**11))**2))
    if lam >0:
        flux=(q_0*np.sin(lam))/(r**2)
    elif lam<0:
        flux=0
    else:
        flux=(q_0*0.0196)/(r**2)
    return flux

def S_32(lam, t, del_0, a, e, phi_ini, m_star, L): #stellar flux for 3-2 resonance
    period=np.sqrt((4*(np.pi**2)*((a*(1.49598*(10**11)))**3))/((6.67259*(10**-11))*m_star))
    E_ini=2*np.arctan(np.sqrt(((1-e)/(1+e)))*np.tan(phi_ini/2))
    M_ini=E_ini-e*np.sin(E_ini)
    d_M=2*np.pi*t*60*60/period
    M=M_ini+d_M
    E=M
    for i in range(8):
        E=M+e*np.sin(E)
    r=a*(1-e*np.cos(E))
    q_0=L/(4*np.pi*((1.49598*(10**11))**2))
    if lam >=(np.pi/2):
        flux=(q_0*0.0523)/(2*np.pi*(r**2))
    elif lam > -(np.pi/2):
        flux=(q_0*np.cos(lam))/(np.pi*(r**2))
    else:
        flux=(q_0*0.0523)/(2*np.pi*(r**2))
    return flux

def A(T): #albedo, depend on temperature as ice cover varies with temperature
    albedo=0.525-0.245*np.tanh((T-268)/5)
    return albedo

def D(omega_p): #diffusitivity: heat transfer assumed mainly by diffusion, summarieses energy transfer processes in the atmosphere
    omega_earth=2*np.pi/(24*60*60)
    diffusitivity_eq=0.5394*((omega_earth/omega_p)**2)
    diffusitivity=0
    return diffusitivity

def C(T,f_o): #heat capacity of surface, influenced by water and ice cover
    C_1=5.25*(10**6)
    C_o=40*C_1
    f_T=(T-273)/10
    if T>=273:
        f_i=0
        Capacity=(1-f_o)*C_1+f_o*((1-f_i)*C_o)
    elif T>=263:
        f_i=1
        C_i=9.2*C_1
        Capacity=(1-f_o)*C_1+f_o*((1-f_i)*C_o+f_i*C_i)
    else:
        f_i=1
        C_i=2*C_1
        Capacity=(1-f_o)*C_1+f_o*((1-f_i)*C_o+f_i*C_i)
    return Capacity

def I(T): #infrared emission from planet
    tau_IR=0.79*((T/273)**3)
    longwave=((5.67037*(10**-8))*T**4)/(1+(3/4)*tau_IR)
    return longwave

#example code to examine behaviour of variables (flux for tidally locked in this case)
#lam=[-90+i*2 for i in range(91)]
#phi_ini=(np.pi/2)
#a=1
#e=0.0167
#flux_11=[S_11(lam[i]*(np.pi/180), 120, a, e, phi_ini, 0.12*1.9885*(10**30), 0.00155*3.84*(10**26)) for i in range(len(lam))]
#plt.plot(lam, flux_11)
#plt.xlabel('latitude(degree)')
#plt.ylabel('S(Wm$^{-2}$)')
#plt.legend()
#plt.show()

