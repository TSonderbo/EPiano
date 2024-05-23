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

    lenSound = fs #length of sound in samples

    rho = 7850 #material Density
    r = 1.524 * 10 **-3
    A = math.pi * r**2 #cross-sectional area
    E = 2*10**11 #young's Modulus
    I = math.pi * r**4 / 4 #inertia

    freq = 41.2

    kappa_1 = math.sqrt(E / rho)
    K = r/2

    L = math.sqrt((1.426 * math.pi * K * kappa_1 / freq) / 8)


    kappa = math.sqrt((E*I)/(rho*A)) #stiffness coefficient

    sigma_0 = 0.001 #freq. independent damping
    sigma_1 = 0.005 #freq. dependent damping

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
    f_alpha = 2.8 #Local geometry of impact coefficient - typically in range [1.5, 3.5]
    f_mu = 0.6
    f_m = 1.1 * 10**-2 #Hammer mass
    f_vIn = 4 #Initial velocity (m/s) typical range [1, 4] - Would be mapped to mpe-midi input
    f_u = -1*10**-4 #Current hammer position
    f_uPrev = f_u - k*f_vIn #Set initial velocity
    f_uNext = 0
    f_ratio = f_m / (rho*A*L)
    f_contact = int(N * 0.4) #math.floor(M_u*0.25)
    x = 0 #Compression
    xPrev = 0 #Previous Compression

    #Damper
    psiNext = np.zeros(N)
    psiPrev = np.zeros(N)
    g = np.zeros(N)
    eta = np.zeros(N)
    etaNext = np.zeros(N)
    etaStar = np.zeros(N)
    etaPrev = np.zeros(N)
    d_kappa = 0 #
    d_K = 1*10**3 #damper stiffness
    d_alpha = 1.3 #damper coefficient
    d_contact = range(int(N*0.8), int(N*0.9)) #contact range
    d_b = 0.00000001 #damper position

    #output
    outLoc = N #ouput location
    output = np.zeros((lenSound,1)) #output vector
    
    f_out = np.zeros((lenSound,1))
    x_out = np.zeros((lenSound,1))

    uPrev = np.copy(u)

    for i in range(0, lenSound):
        
        #Hammer Compression
        x = f_u - u[f_contact]
        
        if(x <= 0): #No contact
            F[f_contact] = 0
            x = 0
        else: #Contact
            #Hunt and cross model
            #F[f_contact] = f_k * x ** f_alpha + f_lambda * x**f_alpha * ((1/k)*(x - xPrev))
            F[f_contact] = f_k * x ** f_alpha * ( 1 + f_mu * ((1/k)*(x - xPrev)))
            
        f_uNext = 2*f_u - f_uPrev - (sum(F)*k**2)
        f_uPrev = f_u
        f_u = f_uNext
        
        xPrev = x
        
        f_out[i] = F[f_contact]
        x_out[i] = x
        
        #Update equation for beam
        for l in range(2, N-1):

            uNext[l] = ((2-6*muSq-((4*sigma_1*k)/hSq))*u[l] \
                + (4*muSq + ((2*sigma_1*k)/hSq))*(u[l+1] + u[l-1]) \
                - muSq*(u[l+2] + u[l-2]) \
                + (-1+sigma_0*k+((4*sigma_1*k)/hSq))*uPrev[l] \
                - ((2*sigma_1*k)/hSq) * (uPrev[l+1] + uPrev[l-1])) \
                + (f_ratio * F[l] * k**2 * 10) \
                /(1+sigma_0*k) 
        
        uStar = np.copy(uNext)
                
        if(i > 0.5*lenSound):
            for l in d_contact:
                
                eta[l] = d_b - u[l]
                etaStar[l] = d_b - uStar[l]
                
                if(psiPrev[l] >= 0):
                    d_kappa = 1
                else:
                    d_kappa = -1
                
                if(eta[l] >= 0):
                    g[l] = d_kappa * math.sqrt((d_K*(d_alpha + 1))/2) * eta[l] ** ((d_alpha -1)/2)
                elif(eta[l] < 0 and etaStar[l] != etaPrev[l]):
                    g[l] = -2*((psiPrev[l])/(etaStar[l]-etaPrev[l]))
                elif(eta[l] < 0 and etaStar[l] == etaPrev[l]):
                    g[l] = 0
                    
                damping = (psiPrev[l]*g[l] + (g[l]**2/4)*uPrev[l])*k**2
                
                uNext[l] = ((2-6*muSq-((4*sigma_1*k)/hSq))*u[l] \
                + (4*muSq + ((2*sigma_1*k)/hSq))*(u[l+1] + u[l-1]) \
                - muSq*(u[l+2] + u[l-2]) + (-1+sigma_0*k+((4*sigma_1*k)/hSq))*uPrev[l] \
                - ((2*sigma_1*k)/hSq) * (uPrev[l+1] + uPrev[l-1]) \
                + damping + (f_ratio * F[l] * k**2)\
                )/ (1+sigma_0*k + (g[l]**2/4)) 
        
        #Update for free boundary based on matrix implementation
        uNext[N-1] = ((2 - 5*muSq - ((4 * sigma_1 * k) / hSq)) * u[N-1] \
        + ((2*sigma_1*k)/(hSq)+2*muSq) * u[N] \
        + (4*muSq + (2*sigma_1*k)/(hSq)) * u[N-2] \
        - muSq*u[N-3]
        + (-1 + sigma_0 * k + (4*sigma_1*k)/hSq)*uPrev[N-1]
        - ((2*sigma_1*k)/(hSq))*(uPrev[N-2] + uPrev[N])) \
        /(1+sigma_0*k)
        
        uNext[N] = ((2 - 2*muSq - ((4 * sigma_1 * k) / hSq)) * u[N] \
        + (4 * muSq + ((4 * sigma_1 * k) / hSq)) * u[N-1] \
        - muSq * 2 * u[N-2] \
        + (-1 + sigma_0 * k + ((4 * sigma_1 * k) / hSq)) * uPrev[N] \
        - ((2 * sigma_1 * k) / hSq) * (2 * uPrev[N-1])) \
        /(1+sigma_0*k)


        for l in d_contact:
            etaNext[l] = uNext[l] - d_b
            psiNext[l] = psiPrev[l] + ((etaNext[l]-etaPrev[l])/2)


        #read output at tip of tine

        output[i] = u[N]

        #update state vectors
        uPrev = np.copy(u)
        u = np.copy(uNext)
        
        etaPrev = np.copy(etaNext)
        psiPrev = np.copy(psiNext)
    
    return output, f_out, x_out
    

@njit
def sim_con():
    
    fs = 48000 * 32 #sample rate
    k = 1/fs #sampling Period
    kSq = k*k

    lenSound = fs #length of sound in samples

    
    rho = 7850 #material Density
    r = 1.524 * 10 **-3
    A = math.pi * r**2 #cross-sectional area
    E = 2*10**11 #young's Modulus
    I = math.pi * r**4 / 4 #inertia


    freq = 41.2

    kappa_1 = math.sqrt(E / rho)
    K = r/2

    L = math.sqrt((1.426 * math.pi * K * kappa_1 / freq) / 8)

    kappa = math.sqrt((E*I)/(rho*A)) #stiffness coefficient

    sigma_0 = 0.001 #freq. independent damping
    sigma_1 = 0.005 #freq. dependent damping

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
    f_k = 1.5 * 10**9 #Hammer stiffness/Spring constant
    f_alpha = 2.8 #Local geometry of impact coefficient - typically in range [1.5, 3.5]
    f_mu = 0.6
    f_m = 1.1 * 10**-2 #Hammer mass
    f_vIn = 4 #Initial velocity (m/s) typical range [1, 4] - Would be mapped to mpe-midi input
    f_u = -1*10**-4 #Current hammer position
    f_uPrev = f_u - k*f_vIn #Set initial velocity
    f_uNext = 0
    f_ratio = f_m / (rho*A*L)
    f_contact = int(N * 0.4) #math.floor(M_u*0.25)
    x = 0 #Compression
    xPrev = 0 #Previous Compression

    #Damper
    d_xi = (0.9 * N)
    d_x = int(d_xi)
    
    K1 = 100
    K3 = 100000
    R = 0.5

    #output
    outLoc = N #ouput location
    output = np.zeros((lenSound,1)) #output vector
    
    f_out = np.zeros((lenSound,1))
    x_out = np.zeros((lenSound,1))

    uPrev = np.copy(u)

    for i in range(0, lenSound):
        
        #Hammer Compression
        x = f_u - u[f_contact]
        
        if(x <= 0): #No contact
            F[f_contact] = 0
            x = 0
        else: #Contact
            #Hunt and cross model
            #F[f_contact] = f_k * x ** f_alpha + f_lambda * x**f_alpha * ((1/k)*(x - xPrev))
            F[f_contact] = f_k * x ** f_alpha * ( 1 + f_mu * ((1/k)*(x - xPrev)))
            
        f_uNext = 2*f_u - f_uPrev - (sum(F)*k**2)
        f_uPrev = f_u
        f_u = f_uNext
        
        xPrev = x
        
        f_out[i] = F[f_contact]
        x_out[i] = x
        
        #Update equation for beam
        for l in range(2, N-1):

            uNext[l] = ((2-6*muSq-((4*sigma_1*k)/hSq))*u[l] \
                + (4*muSq + ((2*sigma_1*k)/hSq))*(u[l+1] + u[l-1]) \
                - muSq*(u[l+2] + u[l-2]) \
                + (-1+sigma_0*k+((4*sigma_1*k)/hSq))*uPrev[l] \
                - ((2*sigma_1*k)/hSq) * (uPrev[l+1] + uPrev[l-1])) \
                + (f_ratio * F[l] * k**2 * 10) \
                /(1+sigma_0*k) 
        
        uStar = np.copy(uNext)
                
        if(0.5 * lenSound < i):

            d_eta = u[d_x]
            d_etaPrev = uPrev[d_x]
            
            wStar = uNext[d_x]
            
            rPlus = K1 / 4 + K3 * (d_eta * d_eta) / 2 + R / (2 * k)
            rMinus = K1 / 4 + K3 * (d_eta * d_eta) / 2 - R / (2 * k)
                
            d_force = ((wStar + (K1 / (2 * rPlus) * d_eta) + (rMinus / (rPlus)*d_etaPrev)) \
		                / (1 / rPlus + (1/h * kSq) / (rho * A * (1 + sigma_0 * k)))) * 0.5
                  
            uNext[d_x] = wStar - (1/h)*((kSq*d_force)/(rho * A * (1+sigma_0*k)))
        
        #Update for free boundary based on matrix implementation
        uNext[N-1] = ((2 - 5*muSq - ((4 * sigma_1 * k) / hSq)) * u[N-1] \
        + ((2*sigma_1*k)/(hSq)+2*muSq) * u[N] \
        + (4*muSq + (2*sigma_1*k)/(hSq)) * u[N-2] \
        - muSq*u[N-3]
        + (-1 + sigma_0 * k + (4*sigma_1*k)/hSq)*uPrev[N-1]
        - ((2*sigma_1*k)/(hSq))*(uPrev[N-2] + uPrev[N])) \
        /(1+sigma_0*k)
        
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

# output, f_out, x_out = sim()
# plt.plot(x_out, f_out)
# plt.show()

# plt.plot(output * 100)
# plt.show()
