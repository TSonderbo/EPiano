import scipy as sp
import numpy as np
import math
import matplotlib.pyplot as plt
from numba import njit
from scipy.fftpack import fft, fftfreq

@njit
def sim():
    
    fs = 48000 * 32 #sample rate
    k = 1/fs #sampling Period

    lenSound = math.floor(fs * 0.5) #length of sound in samples

    L = 0.1564
    rho = 7850 #material Density
    r = 1.524 * 10 **-3
    A = math.pi * r**2 #cross-sectional area
    E = 2*10**11 #young's Modulus
    I = math.pi * r**4 / 4 #inertia

    kappa = math.sqrt((E*I)/(rho*A)) #stiffness coefficient

    sigma_1 = 0.005 #freq. dependent damping
    sigma_0 = 0.0001 #freq. independent damping

    h = math.sqrt((4*sigma_1*k+math.sqrt((4*sigma_1*k)**2+16*kappa**2*k**2))/2) #grid spacing based on stability condition

    N = int(L//h) #grid intervals
    N_frac = L/h

    hSq = h**2

    mu = (kappa*k) / hSq
    muSq = mu**2

    M_u = math.ceil(0.5*N) + 1
    M_w = math.floor(0.5*N) + 1
    
    alpha = N_frac - N

    #state vectors
    uNext = np.zeros(M_u)
    u = np.zeros(M_u)

    wNext = np.zeros(M_w)
    w = np.zeros(M_w)

    #raised cosine
    loc = round(0.5 * M_w); # Center location
    halfWidth = round(N/10); # Half-width of raised cosine
    width = 2 * halfWidth; # Full width
    rcX = range(0,width); # x-locations for raised cosine
    rc = np.array([0.5 - 0.5 * math.cos(2 * math.pi * x / width) for x in rcX]) # raised cosine
    w[loc-halfWidth : loc+halfWidth] = rc # initialise current state

    #output
    outLoc = N #ouput location
    output = np.zeros((lenSound,1)) #output vector

    uPrev = np.copy(u)
    wPrev = np.copy(w)

    for i in range(0, lenSound):
                
        #Update equation for beam
        for l in range(2, N-1):

            uNext[l] = ((2-6*muSq-((4*sigma_1*k)/hSq))*u[l] \
                + (4*muSq + ((2*sigma_1*k)/hSq))*(u[l+1] + u[l-1]) \
                - muSq*(u[l+2] + u[l-2]) \
                + (-1+sigma_0*k+((4*sigma_1*k)/hSq))*uPrev[l] \
                - ((2*sigma_1*k)/hSq) * (uPrev[l+1] + uPrev[l-1])) \
                /(1+sigma_0*k) 
        
        #Update for free boundary based on matrix implementation
        uNext[N-1] = ((2 - 5*muSq - ((4 * sigma_1 * k) / hSq)) * u[N-1] \
        + ((2*sigma_1*k)/(hSq)+2*muSq) * u[N] \
        + (4*muSq + (2*sigma_1*k)/(hSq)) * u[N-2] \
        - muSq*u[N-3]
        + (-1 + sigma_0 * k + (4*sigma_1*k)/hSq)*uPrev[N-1]
        - ((2*sigma_1*k)/(hSq))*(uPrev[N-2] + uPrev[N]))
        
        uNext[N] = ((2 - 2*muSq - ((4 * sigma_1 * k) / hSq)) * u[N] \
        + (4 * muSq + ((4 * sigma_1 * k) / hSq)) * u[N-1] \
        - muSq * 2 * u[N-2] \
        + (-1 + sigma_0 * k + ((4 * sigma_1 * k) / hSq)) * uPrev[N] \
        - ((2 * sigma_1 * k) / hSq) * (2 * uPrev[N-1])) \
        /(1+sigma_0*k)

        #read output at tip of tine
        output[i] = u[outLoc]

        #update state vectors
        uPrev = np.copy(u)
        u = np.copy(uNext)
    
    return output
    

output = sim()

plt.plot(output * 100)
plt.show()
