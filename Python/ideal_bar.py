import scipy as sp
import numpy as np
import math
import matplotlib.pyplot as plt
from numba import njit

@njit
def sim():
    
    fs = 48000 * 16 #sample rate
    k = 1/fs #sampling Period

    lenSound = math.floor(fs * 0.2) #length of sound in samples

    L = 0.1564
    rho = 7850 #material Density
    r = 1.524 * 10 **-3
    A = math.pi * r**2 #cross-sectinal area
    E = 2*10**11 #young's Modulus
    I = math.pi * r**4 / 4 #inertia

    kappa = math.sqrt((E*I)/(rho*A)) #stiffness coefficient

    sigma_1 = 0.005 #freq. dependent damping
    sigma_0 = 0.0001 #freq. independent damping

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
    # loc = round(0.8 * N); # Center location
    # halfWidth = round(N/10); # Half-width of raised cosine
    # width = 2 * halfWidth; # Full width
    # rcX = range(0,width); # x-locations for raised cosine
    # rc = np.array([0.5 - 0.5 * math.cos(2 * math.pi * x / width) for x in rcX]) # raised cosine
    # u[loc-halfWidth : loc+halfWidth] = rc # initialise current state

    #Hammer model
    F = np.zeros(N+1) # Hammer impact force vector - to calculate force across a vector for interpolation
    f_k = 1.5 * 10**11 #Hammer stiffness/Spring constant
    f_alpha = 2.5 #Local geometry of impact coefficient - typically in range [1.5, 3.5]
    f_lambda = (3/2)*f_k #Damping constant
    f_m = 1.1 * 10**-2 #Hammer mass
    f_vIn = 10 #Initial velocity (m/s) typical range [1, 4] - Would be mapped to mpe-midi input
    f_u = -1*10**-4 #Current hammer position
    f_uPrev = f_u - k*f_vIn #Set initial velocity
    f_uNext = 0
    
    f_ratio = f_m / (rho*A*L)
    
    f_contact = math.floor(N*0.8)

    x = 0 #Compression
    xPrev = 0 #Previous Compression

    #output
    outLoc = N #ouput location
    output = np.zeros((lenSound,1)) #output vector
    
    f_out = np.zeros((lenSound,1))
    x_out = np.zeros((lenSound,1))

    uPrev = np.copy(u)

    for i in range(0, lenSound):
        
        #Compression
        x = f_u - u[f_contact]
        
        if(x <= 0): #No contact
            F[f_contact] = 0
        else: #Contact
            #Hunt and cross model
            F[f_contact] = f_k * x ** f_alpha + f_lambda * x**f_alpha * ((1/k)*(x - xPrev))     
            #F[f_contact] = f_k * x ** f_alpha * ( 1 + 0.6* ((1/k)*(x - xPrev)))
            
        f_uNext = 2*f_u - f_uPrev - (F[f_contact]*k**2)/f_m
        f_uPrev = f_u
        f_u = f_uNext
        
        xPrev = x
        
        f_out[i] = F[f_contact]
        x_out[i] = x
        
        #Update equation for beam
        for l in range(2, N-1):
            uNext[l] = ((2-6*muSq-((4*sigma_1*k)/hSq))*u[l] \
                + (4*muSq + ((2*sigma_1*k)/hSq))*(u[l+1] + u[l-1]) \
                - muSq*(u[l+2] + u[l-2]) + (-1+sigma_0*k+((4*sigma_1*k)/hSq))*uPrev[l] \
                - ((2*sigma_1*k)/hSq) * (uPrev[l+1] + uPrev[l-1])) \
                /(1+sigma_0*k) + (f_ratio * F[f_contact] * k**2)
        
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

        output[i] = u[N]

        #update state vectors
        uPrev = np.copy(u)
        u = np.copy(uNext)
    
    return output, f_out, x_out
    

output, f_out, x_out = sim()
plt.plot(x_out, f_out)
plt.show()

plt.plot(output)
plt.show()