import scipy as sp
import numpy as np
import math
import matplotlib.pyplot as plt
from numba import njit

@njit
def sim():
    
    fs = 48000 * 4 #sample rate
    k = 1/fs #sampling Period

    lenSound = 4 * fs #length of sound in samples

    L = 0.1564
    rho = 7850 #material Density
    r = 0.001905 / 2
    A = math.pi * r**2 #cross-sectinal area
    E = 2*10**11 #young's Modulus
    I = math.pi * r**4 / 4 #inertia

    kappa = math.sqrt((E*I)/(rho*A)) #stiffness coefficient

    sigma_1 = 0.05 #freq. dependent damping
    sigma_0 = 1.0 #freq. independent damping

    h = math.sqrt((4*sigma_1*k+math.sqrt((4*sigma_1*k)**2+16*kappa**2*k**2))/2) #grid spacing based on stability condition

    N = int(L//h) #grid intervals

    h = L/N 
    hSq = h**2

    mu = (kappa*k) / hSq
    muSq = mu**2

    #state vectors
    uNext = np.zeros(N+1)
    u = np.zeros(N+1)
    vg1Prev = 0

    # raised cosine
    loc = round(0.5 * N); # Center location
    halfWidth = round(N/10); # Half-width of raised cosine
    width = 2 * halfWidth; # Full width
    rcX = range(0,width); # x-locations for raised cosine
    rc = np.array([0.5 - 0.5 * math.cos(2 * math.pi * x / width) for x in rcX]) # raised cosine
    u[loc-halfWidth : loc+halfWidth] = rc # initialise current state

    #output
    outLoc = N #ouput location
    output = np.zeros((lenSound,1)) #output vector

    uPrev = np.copy(u)

    for i in range(0, lenSound):
        for l in range(2, N-1):
            uNext[l] = ((2-6*muSq-((4*sigma_1*k)/hSq))*u[l] \
                + (4*muSq + ((2*sigma_1*k)/hSq))*(u[l+1] + u[l-1]) \
                - muSq*(u[l+2] + u[l-2]) + (-1+sigma_0*k+((4*sigma_1*k)/hSq))*uPrev[l] \
                - ((2*sigma_1*k)/hSq) * (uPrev[l+1] + uPrev[l-1])) \
                /(1+sigma_0*k)   

        #boundary conditions - clamped/free
        #virtual grid points for free boundary
        vg1 = 2*u[N] - u[N-1] # N+1
        vg2 = u[N-2]-2*u[N-1]+2*vg1 # N+2

        uNext[N-1] = ((2-6*muSq-(4*sigma_1*k)/hSq)*u[N-1] \
            + (4*muSq + (2*sigma_1*k)/hSq)*(u[N] + u[N-2]) \
            - muSq*(vg1 + u[N-3]) + (-1+sigma_0*k+(4*sigma_1*k)/hSq)*uPrev[N-1] \
            - ((2*sigma_1*k)/hSq) * (uPrev[N] + uPrev[N-2]))/(1+sigma_0*k)   

        uNext[N] = ((2-6*muSq-(4*sigma_1*k)/hSq)*u[N] \
            + (4*muSq + (2*sigma_1*k)/hSq)*(vg1 + u[N-1]) \
            - muSq*(vg2 + u[N-2]) + (-1+sigma_0*k+(4*sigma_1*k)/hSq)*uPrev[N] \
            - ((2*sigma_1*k)/hSq) * (vg1Prev + uPrev[N-1]))/(1+sigma_0*k)   

        vg1Prev = vg1

        #read output at tip of tine

        output[i] = u[N]

        #update state vectors
        uPrev = np.copy(u)
        u = np.copy(uNext)
    
    return output
    

output = sim()
plt.plot(output)
plt.show()