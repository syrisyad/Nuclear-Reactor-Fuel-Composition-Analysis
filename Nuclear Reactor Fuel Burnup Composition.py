# Parameter

Fraksi = 0.5
T = 1000   # Batas waktu (tahun)
dT = 1
Flux = (3.11E+13)/3.17098e-8 # (Neutron/cm tahun)

# Penampang lintang penyerapan (barn)

CSa_u238 = 12.0E-24
CSa_pu239 = 972.0E-24
CSa_pu240 = 317.0E-24
CSa_pu241 = 1264.0E-24
CSa_pu242 = 17.0E-24

# Penampang lintang penangkapan (barn)

CSc_u238 = 2.0E-24
CSc_pu239 = 274.0E-24
CSc_pu240 = 264.0E-24
CSc_pu241 = 326.0E-24

# Konstanta peluruhan (1/tahun)

Lambda_u239 = 4.92E-4/3.17098e-8
Lambda_np239 = 34.0E-9/3.17098e-8
Lambda_pu241 = 1.53E-9/3.17098e-8

import numpy as np

L = int(T/dT)

# Inisiasi array

waktu = np.zeros(L)
t = 0.0
N_u238  = np.zeros(L)
N_u239  = np.zeros(L)
N_np239 = np.zeros(L)
N_pu239 = np.zeros(L)
N_pu240 = np.zeros(L)
N_pu241 = np.zeros(L)
N_pu242 = np.zeros(L)

# Banyak neutron mula

N_u238[0] = 5.0E+20
N_u239[0] = 0
N_np239[0] = 0
N_pu239[0] = 0
N_pu240[0] = 0
N_pu241[0] = 0
N_pu242[0] = 0

# Persamaan semi-implisit

for i in range (1,L):
    waktu[i] = 0+t
    N_u238[i] = N_u238[i-1]*(1-(CSa_u238*Flux*dT*(1-Fraksi)))/(1+(CSa_u238*Flux*Fraksi*dT))
    N_u239[i] = (CSc_u238*Flux*dT*((N_u238[i]*Fraksi)+(N_u238[i-1]*(1-Fraksi)))+N_u239[i-1]*(1-(Lambda_u239*(1-Fraksi)*dT)))/(1+(Lambda_u239*Fraksi*dT))
    N_np239[i] = (Lambda_u239*dT*((N_u239[i]*Fraksi)+(N_u239[i-1]*(1-Fraksi)))+(N_np239[i-1]*(1+(Lambda_np239*(1-Fraksi)*dT))))/(1+Lambda_np239*Fraksi*dT)
    N_pu239[i] = (Lambda_np239*dT*(N_np239[i]*Fraksi+N_np239[i-1]*(1-Fraksi))+N_pu239[i-1]*(1-CSa_pu239*Flux*(1-Fraksi)*dT))/(1+CSa_pu239*Flux*Fraksi*dT)
    N_pu240[i] = (CSc_pu239*Flux*dT*(N_pu239[i]*Fraksi+N_pu239[i-1]*(1-Fraksi))+N_pu240[i-1]*(1-CSa_pu240*Flux*dT*(1-Fraksi)))/(1+CSa_pu240*Flux*Fraksi*dT)
    N_pu241[i] = (CSc_pu240*Flux*dT*(N_pu240[i]*Fraksi+N_pu240[i-1])+N_pu241[i-1]*(1-Fraksi*dT*(CSa_pu241*Flux+Lambda_pu241)))/(1+Fraksi*dT*(CSa_pu241*Flux+Lambda_pu241))
    N_pu242[i] = (CSc_pu241*Flux*dT*(N_pu241[i]*Fraksi+N_pu241[i-1])+N_pu242[i-1]*(1-CSa_pu242*Flux*dT*(1-Fraksi)))/(1+CSa_pu242*Flux*dT*Fraksi)
    t = t + dT

import matplotlib.pyplot as plt

#Plot semilog

plt.semilogy(waktu,N_u238, label = 'U-238')
plt.semilogy(waktu,N_pu239, label ='Pu-239')
plt.semilogy(waktu,N_pu242, label ='Pu-242')

plt.legend()
plt.title('Grafik Burnup UO2, Fraksi 0.5')
plt.xlabel('Waktu(tahun)')
plt.ylabel('Neutron')

plt.grid(True, which="both")

plt.show()
