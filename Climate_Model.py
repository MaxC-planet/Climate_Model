import numpy as np
import Get_SADCI as get
import matplotlib.pyplot as plt
#import temporal_fractional_hab as tfh (not used in project)

# set length of time simulated, time step and latitude band thickness, and also host star mass and luminosity
L_star=0.00155*3.84*(10**26)
m_star=0.12*1.9885*(10**30)
t_tot=150*60*60*24*365 #year
dt=0.02*60*60*24 #day
band_size=9 # in degrees, choose values such that band no is already integer!

tstep_no=int(t_tot/dt)+1
band_no=int(180/band_size)+1 # 2 bands surronding poles

band_size_rad=band_size*(np.pi/180) # band size in radians
la=np.array([(-90+band_size*i)*(np.pi/180) for i in range(band_no)])
lam=np.zeros(band_no)
for i in range(band_no):
    lam[i]=la[-i-1]

#initial variables and units:
#del_0: obliquity, in rad
#a: semi-major axis, in AU
#period of rotation, in hrs
    
def central_T11(initial_T, a, e, phi_ini, f_o, period, m_star, L): #all central
    T=np.zeros((band_no,tstep_no))
    omega_p=(2*np.pi)/(period*60*60)
    T[:,0]=initial_T #initial conditions
    for j in range (0,tstep_no-1):
        T[0,j+1]=(dt/get.C(T[0,j],f_o))*(get.S_11(lam[0],j*dt,a,e,phi_ini,m_star,L)*(1-get.A(T[0,j]))+
        get.D(omega_p)*(T[2,j]-T[0,j])/(2*(band_size_rad**2))-
        get.I(T[0,j]))+T[0,j]
        T[-1,j+1]=(dt/get.C(T[-1,j],f_o))*(get.S_11(lam[-1],j*dt,a,e,phi_ini,m_star,L)*(1-get.A(T[-1,j]))+
        get.D(omega_p)*(T[-3,j]-T[-1,j])/(2*(band_size_rad**2))-
        get.I(T[-1,j]))+T[-1,j]
        for i in range (1,band_no-1):
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

def central_T32(initial_T, del_0, a, e, phi_ini, f_o, period, m_star, L): #all central
    T=np.zeros((band_no,tstep_no), dtype='float64')
    omega_p=(2*np.pi)/(period*60*60)
    T[:,0]=initial_T #initial conditions
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

period_p=11.186 #day, also 7.457
#T1=central_T32(np.array([300 for i in range(band_no)]), 23.436*(np.pi/180), 1, 0.0167, (np.pi/2), 0.7, 24, m_star, L_star) 
#T2=central_T32(np.array([300 for i in range(band_no)]), 0, 0.0485, 0.35, 0, 0.4, 24*7.457, m_star, L_star)
#T3=central_T32(np.array([300 for i in range(band_no)]), 0, 0.0485, 0.35, 0, 0.7, 24*7.457, m_star, L_star)
#T4=central_T32(np.array([300 for i in range(band_no)]), 0, 0.0485, 0.35, 0, 0.85, 24*7.457, m_star, L_star)
#T5=central_T32(np.array([300 for i in range(band_no)]), 0, 0.0485, 0.35, 0, 1, 24*7.457, m_star, L_star)
lat=[90,75,60,45,30,15,0,-15,-30,-45,-60,-75,-90]
band_1day=[5,6,7.5,9,10]
band_2day=[7.5,9,10]
band_3day=[9,10]
T_1day=[293.5272317918473,292.01874070742264,291.52904482213427,291.0314190561704,290.1675886536959]
T_2day=[292.2537964471909,291.2173659532141,290.377035586764]
T_3day=[291.5699168106933,290.5602763397645]
#plt.plot(band_2day, T_2day, label='2d')
#plt.plot(band_3day, T_3day, label='3d')
#plt.xlabel('latitude(degrees)')
#plt.ylabel('global mean temp(K)')
#plt.legend()

#t=[i for i in range(len(T1[0,:])-500*365)]
#lat=[90,45,0,-45,-90]
#for i in range(len(lat)):
    #plt.plot(t, T1[i*5,500*365:],label=lat[i])
#plt.xticks([i*100*365*10 for i in range (10)], [i*10+5 for i in range (10)])
#plt.ylabel('T(K)')
#plt.xlabel('t(yr)')
#plt.legend()
#print(la)
#print(lam)
    
#mean_T1=np.mean(T1[:,500*365:], axis=1)
#mean_T2=np.mean(T2[:,500*365:], axis=1)
#mean_T3=np.mean(T3[:,500*365:], axis=1)
#mean_T4=np.mean(T4[:,500*365:], axis=1)
#mean_T5=np.mean(T5[:,500*365:], axis=1)
parametrised_pro=[288-30.2*((3*((np.sin(lam[i]))**2)-1)/2) for i in range(len(lam))]
#print(GWMT.GT(lam, parametrised_pro, band_size_rad))
#print(GWMT.GT(lam, mean_T, band_size_rad))
#plt.scatter(lam*(180/np.pi), mean_T1, label='model')
#plt.plot(lam*(180/np.pi), parametrised_pro, label='parametrized')
#plt.plot(lam*(180/np.pi), mean_T1, label='0.1')
#plt.plot(lam*(180/np.pi), mean_T2, label='0.4')
#plt.plot(lam*(180/np.pi), mean_T3, label='0.7')
#plt.plot(lam*(180/np.pi), mean_T4, label='0.85')
#plt.plot(lam*(180/np.pi), mean_T5, label='1')
lab=np.array([0.06,0.13,0.2,0.27,0.35])
spin=(np.sqrt((4*(np.pi**2)*((lab*(1.49598*(10**11)))**3))/((6.67259*(10**-11))*m_star)))/(60*60)
for i in range(5):
    plt.plot(lam*(180/np.pi), np.mean(central_T32(np.array([300 for i in range(band_no)]), 0, 0.0485, lab[i], 0, 0.7, 24*7.457, m_star, L_star)[:,2500*365:], axis=1), label=lab[i])
#plt.title('Temperature profile for different ocean fractions')
#plt.yticks([255+i*5 for i in range(11)])
plt.ylabel('T(K)')
plt.xlabel('\u03bb(degree)')
plt.legend()
#print(np.sum((mean_T-np.array(parametrised_pro))**2))
#print(mean_T)

#plt.imshow(T1, interpolation='bicubic', aspect='auto', cmap='hot')
#plt.xticks([i*50*365 for i in range(11)], [i*50 for i in range(11)])
#plt.yticks([i*2 for i in range(11)], [90-i*18 for i in range(11)])
#plt.xlabel('time(year)')
#plt.ylabel('latitude(degree)')
#plt.colorbar()
#plt.show()
#plt.legend()

#sub.run(['echo','getting array'])
ran_e_11=7
ran_e_32=30
ran_fo=20
#tem_h_map=np.zeros((ran_e_11, ran_fo+1)) #change between 32 and 11
#for i in range (ran_e_11):
    #tem_h_map[i,0]=tfh.tem_h(central_T11(np.array([300 for i in range(band_no)]), 0.0485, 0+i*0.01, 0, 0.01, 24*11.186, m_star, L_star)[:,1000*365:])
#for i in range(ran_e_11):
    #for j in range (ran_fo): #change between 32 and 11
        #tem_h_map[i,j+1]=tfh.tem_h(central_T11(np.array([300 for i in range(band_no)]), 0.0485, 0+i*0.01, 0, 0.05+j*0.05, 24*11.186, m_star, L_star)[:,1000*365:]) # change between 32 and 11
#print('===================')
#plt.imshow(tem_h_map, aspect='auto')
#plt.xticks([i*0.01 for i in range(ran_e_11)])
#plt.yticks([0,19,39,59,79,99],[0.01,0.2,0.4,0.6,0.8,1])
#plt.xlabel('eccentricity')
#plt.ylabel('ocean fraction')
#print(tfh.tem_h(T1))
#print(tfh.tem_h(T2))
#sub.run(['echo','getting array'])

#np.savetxt('1:1.out',tem_h_map)
plt.savefig('corrected_pcba_3_e2_tdiff')
#plt.savefig('test_plot.png')
#optimal: 9 degree & 3day
#optimal: 7 degree & 2day  1000yr-5min
#optimal: 5/6 degree & 1day  20yr-less than 1min, 100yr-1min, 1000yr-9min
#optimal: 4 degree & 0.5day  20yr-1min, 100yr-3min, 1000yr-22min
#optimal: 3 degree & 0.1day  20yr-4min, 100yr-16min, 1000yr-2 hr 8min
#optimal: 2 degree & 0.1day  20yr-5min, 100yr-24min, 1000yr-3hr 21min
#higher band size or smaller time step cause smaller spike
#time step of less then few days may not physcially make sense (although less
  #can be used to reduce spike)

#long term cycle of temp fluctuations
#small variations of annual low (and high) temp for time step 2 and 3 day

#optimising model: similar relationship to spike size
#current: 1day, 9 degree

#plotting of global mean temp confirms relationship of dt and dlambda to fit
#dlambda bigger effect than dt on global temp
#will use 1day 9 degree: difference in global mean temp around 3.4K
#bicubic interpolation supresses logn term cycle

#flux of equator mainly affected by obliquity (not affected by vernal equinox time)
#higher obliquity higher temp: can increase habitable a with higher obliquity?
#can do a&e maps for different obliquities

#low ocean fraction can cause snowball state (<0.12), mainly affects position of profile (relative temp between bands largely unaffected-linked to assumption?)
#luminosity slightly below Sun can cause snowball state (<0.99L_Sun!)
#time step for proxima b need small as 0.02d to eliminate errors for 3:2, 0.01d for 1:1 (more like avoid overflow error but no comp errors otherwise)
