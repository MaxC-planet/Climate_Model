import numpy as np
import Get_SADCI as get #module used to model variables in the governing equation, see Get_SADCI.py
import matplotlib.pyplot as plt

#set length of time simulated (t_tot), time step (dt) and latitude band size (band_size)
#also host star mass (m_star) and luminosity (L_star)
#star characteristics taken from observed data on Earth and Proxima Centauri b
L_star=0.00155*3.84*(10**26) #in W
m_star=0.12*1.9885*(10**30) #in kg
t_tot=150*60*60*24*365 #in yr
dt=0.02*60*60*24 #in days
band_size=9 # in degrees, choose values such that total number of bands is already integer

tstep_no=int(t_tot/dt)+1
band_no=int(180/band_size)+1 # 2 bands surronding poles

band_size_rad=band_size*(np.pi/180) # band size in radians
la=np.array([(-90+band_size*i)*(np.pi/180) for i in range(band_no)])
lam=np.zeros(band_no)
for i in range(band_no):
    lam[i]=la[-i-1]

#initial variables and units:
#initial_T: inital temperature, set to 300 K 
#a: semi-major axis, in AU
#e: eccentricity
#phi_ini: inital phase in orbit, in radians, represent inital position in orbit (found to have little effect on results)
#f_o: ocean fraction
#period: period of rotation, in hr
#del_0: obliquity, in radians (only relevant for 3-2 resonance)

#Research found that Proxima Centauri b can have 1-1 (tidally locked) or 3-2 spin-orbit resonance
#11 for tidally locked, 32 for 3-2 resonance
#central difference method for integration
    
def central_T11(initial_T, a, e, phi_ini, f_o, period, m_star, L): 
    T=np.zeros((band_no,tstep_no))
    omega_p=(2*np.pi)/(period*60*60)
    T[:,0]=initial_T #initial conditions
    for j in range (0,tstep_no-1): #evolve temperature distribution 
        T[0,j+1]=(dt/get.C(T[0,j],f_o))*(get.S_11(lam[0],j*dt,a,e,phi_ini,m_star,L)*(1-get.A(T[0,j]))+
        get.D(omega_p)*(T[2,j]-T[0,j])/(2*(band_size_rad**2))-
        get.I(T[0,j]))+T[0,j]
        T[-1,j+1]=(dt/get.C(T[-1,j],f_o))*(get.S_11(lam[-1],j*dt,a,e,phi_ini,m_star,L)*(1-get.A(T[-1,j]))+
        get.D(omega_p)*(T[-3,j]-T[-1,j])/(2*(band_size_rad**2))-
        get.I(T[-1,j]))+T[-1,j]
        for i in range (1,band_no-1): #temperature distribution in one time step
            T[i,j+1]=(dt/get.C(T[i,j],f_o))*(get.S_11(lam[i],j*dt,a,e,phi_ini,m_star,L)*(1-get.A(T[i,j]))-
            get.D(omega_p)*np.tan(lam[i])*(T[i-1,j]-T[i+1,j])/(2*band_size_rad)+
            get.D(omega_p)*(T[i+1,j]-2*T[i,j]+T[i-1,j])/band_size_rad**2-
            get.I(T[i,j]))+T[i,j]
    T[0,tstep_no-1]=(dt/get.C(T[0,tstep_no-2],f_o))*(get.S_11(lam[0],(tstep_no-2)*dt,a,e,phi_ini,m_star,L)*(1-get.A(T[0,tstep_no-2]))+
    get.D(omega_p)*(T[2,tstep_no-2]-T[0,tstep_no-2])/(2*(band_size_rad**2))-
    get.I(T[0,tstep_no-2]))+T[0,tstep_no-2]
    T[-1,tstep_no-1]=(dt/get.C(T[-1,tstep_no-2],f_o))*(get.S_11(lam[-1],(tstep_no-2)*dt,a,e,phi_ini,m_star,L)*(1-get.A(T[-1,tstep_no-2]))+
    get.D(omega_p)*(T[-3,tstep_no-2]-T[-1,tstep_no-2])/(2*(band_size_rad**2))-
    get.I(T[-1,tstep_no-2]))+T[-1,tstep_no-2]
    return T

#3-2 model is also adapted to model Earth
def central_T32(initial_T, del_0, a, e, phi_ini, f_o, period, m_star, L): 
    T=np.zeros((band_no,tstep_no), dtype='float64')
    omega_p=(2*np.pi)/(period*60*60)
    T[:,0]=initial_T 
    for j in range (tstep_no-1):
        T[0,j+1]=(dt/get.C(T[0,j],f_o))*(get.S_32(lam[0],j*dt,del_0,a,e,phi_ini,m_star,L)*(1-get.A(T[0,j]))+
        get.D(omega_p)*(T[2,j]-T[0,j])/(2*(band_size_rad**2))-
        get.I(T[0,j]))+T[0,j]
        T[-1,j+1]=(dt/get.C(T[-1,j],f_o))*(get.S_32(lam[-1],j*dt,del_0,a,e,phi_ini,m_star,L)*(1-get.A(T[-1,j]))+
        get.D(omega_p)*(T[-3,j]-T[-1,j])/(2*(band_size_rad**2))-
        get.I(T[-1,j]))+T[-1,j]
        for i in range (1,band_no-1):
            T[i,j+1]=(dt/get.C(T[i,j],f_o))*(get.S_32(lam[i],j*dt,del_0,a,e,phi_ini,m_star,L)*(1-get.A(T[i,j]))-
            get.D(omega_p)*np.tan(lam[i])*(T[i-1,j]-T[i+1,j])/(2*band_size_rad)+
            get.D(omega_p)*(T[i+1,j]-2*T[i,j]+T[i-1,j])/band_size_rad**2-
            get.I(T[i,j]))+T[i,j]
    T[0,tstep_no-1]=(dt/get.C(T[0,tstep_no-2],f_o))*(get.S_32(lam[0],(tstep_no-2)*dt,del_0,a,e,phi_ini,m_star,L)*(1-get.A(T[0,tstep_no-2]))+
    get.D(omega_p)*(T[2,tstep_no-2]-T[0,tstep_no-2])/(2*(band_size_rad**2))-
    get.I(T[0,tstep_no-2]))+T[0,tstep_no-2]
    T[-1,tstep_no-1]=(dt/get.C(T[-1,tstep_no-2],f_o))*(get.S_32(lam[-1],(tstep_no-2)*dt,del_0,a,e,phi_ini,m_star,L)*(1-get.A(T[-1,tstep_no-2]))+
    get.D(omega_p)*(T[-3,tstep_no-2]-T[-1,tstep_no-2])/(2*(band_size_rad**2))-
    get.I(T[-1,tstep_no-2]))+T[-1,tstep_no-2]
    return T

#(unrefined) code for a temperature distribtuion plot, below example is for different eccentricity values
#mean_T1=np.mean(T1[:,500*365:], axis=1)
#mean_T2=np.mean(T2[:,500*365:], axis=1)
#mean_T3=np.mean(T3[:,500*365:], axis=1)
#mean_T4=np.mean(T4[:,500*365:], axis=1)
#mean_T5=np.mean(T5[:,500*365:], axis=1)
#parametrised_pro=[288-30.2*((3*((np.sin(lam[i]))**2)-1)/2) for i in range(len(lam))]
#plt.plot(lam*(180/np.pi), parametrised_pro, label='parametrized')
#plt.plot(lam*(180/np.pi), mean_T1, label='0.1')
#plt.plot(lam*(180/np.pi), mean_T2, label='0.4')
#plt.plot(lam*(180/np.pi), mean_T3, label='0.7')
#plt.plot(lam*(180/np.pi), mean_T4, label='0.85')
#plt.plot(lam*(180/np.pi), mean_T5, label='1')
#lab=np.array([0.06,0.13,0.2,0.27,0.35])
#for i in range(5):
    #plt.plot(lam*(180/np.pi), np.mean(central_T32(np.array([300 for i in range(band_no)]), 0, 0.0485, lab[i], 0, 0.7, 24*7.457, m_star, L_star)[:,2500*365:], axis=1), label=lab[i])
#plt.title('Temperature profile for different ocean fractions')
#plt.yticks([255+i*5 for i in range(11)])
#plt.ylabel('T(K)')
#plt.xlabel('\u03bb(degree)')
#plt.legend()


#example code for plotting a heat map
#plt.imshow(T1, interpolation='bicubic', aspect='auto', cmap='hot')
#plt.xticks([i*50*365 for i in range(11)], [i*50 for i in range(11)])
#plt.yticks([i*2 for i in range(11)], [90-i*18 for i in range(11)])
#plt.xlabel('time(year)')
#plt.ylabel('latitude(degree)')
#plt.colorbar()
#plt.show()
#plt.legend()

