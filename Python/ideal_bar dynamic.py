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
    kSq = k**2
    
    lenSound = math.floor(fs * 0.5) #length of sound in samples

    L = 0.0864
    rho = 7850 #material Density
    r = 1.524 * 10 **-3
    A = math.pi * r**2 #cross-sectional area
    E = 2*10**11 #young's Modulus
    I = math.pi * r**4 / 4 #inertia

    kappa = math.sqrt((E*I)/(rho*A)) #stiffness coefficient

    sigma_0 = 0.001 #freq. independent damping
    sigma_1 = 0.005 #freq. dependent damping

    h = math.sqrt((4*sigma_1*k+math.sqrt((4*sigma_1*k)**2+16*kappa**2*k**2))/2) #grid spacing based on stability condition

    N = int(L//h) #grid intervals
    N_frac = L/h
    alpha = N_frac - N

    hSq = h**2

    mu = (kappa*k) / hSq
    muSq = mu**2

    M_u = math.ceil(0.5*N)
    M_w = math.floor(0.5*N)

    Iterm = (alpha - 1.0) / (alpha + 1.0)
    itermSq = Iterm * Iterm
    J = np.zeros(4)
    J[0] = itermSq - 4.0 * Iterm + 6.0
    J[1] = Iterm - 4.0
    J[2] = -itermSq + 4.0 * Iterm + 1.0
    J[3] = -Iterm
    
    #state vectors
    uNext = np.zeros(M_u + 1)
    u = np.zeros(M_u + 1)
    uPrev = np.zeros(M_u + 1)
   
    wNext = np.zeros(M_w + 1)
    w = np.zeros(M_w + 1)
    wPrev = np.zeros(M_w + 1)
    
    #Hammer model
    F = np.zeros(M_u+1) # Hammer impact force vector - to calculate force across a vector for interpolation
    f_k = 1.5 * 10**11 #Hammer stiffness/Spring constant
    f_alpha = 2.8 #Local geometry of impact coefficient - typically in range [1.5, 3.5]
    f_mu = 0.6
    f_m = 1.1 * 10**-2 #Hammer mass
    f_vIn = 4 #Initial velocity (m/s) typical range [1, 4] - Would be mapped to mpe-midi input
    f_u = -1*10**-4 #Current hammer position
    f_uPrev = f_u - k*f_vIn #Set initial velocity
    f_uNext = 0
    f_ratio = f_m / (rho*A*L)
    f_contact = 4 #math.floor(M_u*0.25)
    x = 0 #Compression
    xPrev = 0 #Previous Compression
    
    f_out = np.zeros((lenSound,1))
    x_out = np.zeros((lenSound,1))


    # # raised cosine
    # loc = round(0.5 * M_u); # Center location
    # halfWidth = round(M_u/10); # Half-width of raised cosine
    # width = 2 * halfWidth; # Full width
    # rcX = range(0,width); # x-locations for raised cosine
    # rc = np.array([0.5 - 0.5 * math.cos(2 * math.pi * x / width) for x in rcX]) # raised cosine
    # u[loc-halfWidth : loc+halfWidth] = rc # initialise current state
    
    d_xi = (0.8 * M_w)
    d_x = int(d_xi)
    
    K1 = 100
    K3 = 100000
    R = 0.5

    uPrev = u.copy()
    
    #output
    outLoc = int(M_w * 0.5) #ouput location
    output = np.zeros((lenSound,1)) #output vector
    force_output = np.zeros((lenSound,1)) #output vector

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
        for l in range(2, M_w-1):

            uNext[l] = ((2-6*muSq-((4*sigma_1*k)/hSq))*u[l] \
                + (4*muSq + ((2*sigma_1*k)/hSq))*(u[l+1] + u[l-1]) \
                - muSq*(u[l+2] + u[l-2]) \
                + (-1+sigma_0*k+((4*sigma_1*k)/hSq))*uPrev[l] \
                - ((2*sigma_1*k)/hSq) * (uPrev[l+1] + uPrev[l-1])) \
                + (f_ratio * F[l] * k**2 * 10) \
                /(1+sigma_0*k) 
                
            wNext[l] = ((2-6*muSq-((4*sigma_1*k)/hSq))*w[l] \
                + (4*muSq + ((2*sigma_1*k)/hSq))*(w[l+1] + w[l-1]) \
                - muSq*(w[l+2] + w[l-2]) \
                + (-1+sigma_0*k+((4*sigma_1*k)/hSq))*wPrev[l] \
                - ((2*sigma_1*k)/hSq) * (wPrev[l+1] + wPrev[l-1])) \
                /(1+sigma_0*k) 
        
        if(M_u > M_w):
            l = M_u - 2
            
            uNext[l] = ((2-6*muSq-((4*sigma_1*k)/hSq))*u[l] \
                + (4*muSq + ((2*sigma_1*k)/hSq))*(u[l+1] + u[l-1]) \
                - muSq*(u[l+2] + u[l-2]) \
                + (-1+sigma_0*k+((4*sigma_1*k)/hSq))*uPrev[l] \
                - ((2*sigma_1*k)/hSq) * (uPrev[l+1] + uPrev[l-1])) \
                /(1+sigma_0*k) 

        uNext[M_u - 1] = ((2 - 6*muSq - (4*sigma_1*k)/hSq)*u[M_u - 1] \
                          + (-J[1]*muSq + (2*sigma_1*k)/hSq) * u[M_u] \
                          + (4*muSq + (2*sigma_1*k)/hSq) * u[M_u - 2] \
                          - muSq*u[M_u - 3] \
                          - muSq*w[M_w] \
                          - J[3]*muSq*w[M_w - 1] \
                          + (-1+sigma_0*k + (4*sigma_1*k)/hSq) * uPrev[M_u - 1] \
                          - ((2*sigma_1*k)/hSq) * (uPrev[M_u - 2] + uPrev[M_u]) \
                          ) / (1+sigma_0*k) 
        
        uNext[M_u] = ((2 - J[0] * muSq + (Iterm - 2)*((2*sigma_1*k)/hSq)) * u[M_u] \
                      + (4*muSq + ((2*sigma_1*k)/hSq)) * u[M_u - 1] \
                      - muSq * u[M_u - 2] \
                      + (((2*sigma_1*k)/hSq) - J[1]*muSq)*w[M_w] \
                      + (J[3]*((2*sigma_1*k)/hSq) - J[2]*muSq)*w[M_w - 1] \
                      - (J[3]*muSq)*w[M_w - 2] \
                      + (-1 + sigma_0*k - (Iterm - 2)*((2*sigma_1*k)/hSq))*uPrev[M_u] \
                      - ((2*sigma_1*k)/hSq)*(uPrev[M_u - 1] + wPrev[M_w]) \
                      - (J[3] * ((2*sigma_1*k)/hSq)) * wPrev[M_w - 1]\
                      ) / (1+sigma_0*k) 
        
        wNext[M_w - 1] = ((2 - 6*muSq - (4*sigma_1*k)/hSq)*w[M_w - 1] \
                          + (-J[1]*muSq + (2*sigma_1*k)/hSq) * w[M_w] \
                          + (4*muSq + (2*sigma_1*k)/hSq) * w[M_w - 2] \
                          - muSq*w[M_w - 3] \
                          - muSq*u[M_u] \
                          - J[3]*muSq*u[M_u - 1] \
                          + (-1+sigma_0*k + (4*sigma_1*k)/hSq) * wPrev[M_w - 1] \
                          - ((2*sigma_1*k)/hSq) * (wPrev[M_w - 2] + wPrev[M_w]) \
                          ) / (1+sigma_0*k) 
        
        wNext[M_w] = ((2 - J[0] * muSq + (Iterm - 2)*((2*sigma_1*k)/hSq)) * w[M_w] \
                      + (4*muSq + ((2*sigma_1*k)/hSq)) * w[M_w - 1] \
                      - muSq * w[M_w - 2] \
                      + (((2*sigma_1*k)/hSq) - J[1]*muSq)*u[M_u] \
                      + (J[3]*((2*sigma_1*k)/hSq) - J[2]*muSq)*u[M_u - 1] \
                      - (J[3]*muSq)*u[M_u - 2] \
                      + (-1 + sigma_0*k - (Iterm - 2)*((2*sigma_1*k)/hSq))*wPrev[M_w] \
                      - ((2*sigma_1*k)/hSq)*(wPrev[M_w - 1] + uPrev[M_u]) \
                      - (J[3] * ((2*sigma_1*k)/hSq)) * uPrev[M_u - 1]\
                      ) / (1+sigma_0*k) 
        
        #Update for free boundary based on matrix implementation
        wNext[1] = ((2 - 5*muSq - ((4 * sigma_1 * k) / hSq)) * w[1] \
        + ((2*sigma_1*k)/(hSq)+2*muSq) * w[0] \
        + (4*muSq + (2*sigma_1*k)/(hSq)) * w[2] \
        - muSq*w[3] \
        + (-1 + sigma_0 * k + (4*sigma_1*k)/hSq)*wPrev[1] \
        - ((2*sigma_1*k)/(hSq))*(wPrev[2] + wPrev[0]))\
        /(1+sigma_0*k) 
        
        wNext[0] = ((2 - 2*muSq - ((4 * sigma_1 * k) / hSq)) * w[0] \
        + (4 * muSq + ((4 * sigma_1 * k) / hSq)) * w[1] \
        - muSq * 2 * w[2] \
        + (-1 + sigma_0 * k + ((4 * sigma_1 * k) / hSq)) * wPrev[0] \
        - ((2 * sigma_1 * k) / hSq) * (2 * wPrev[1])) \
        /(1+sigma_0*k)

        #apply damper
        if(0.5 * lenSound < i):
            # for l in range(0, 4):
            #     d_eta = w[d_x - 1 + l] * d_I[l] #Distance from 0
            #     d_etaPrev = wPrev[d_x - 1 + l] * d_I[l]
                
            #     wStar = wNext[d_x - 1 + l] * d_I[l]
                
            #     rPlus = K1 / 4 + K3 * (d_eta * d_eta) / 2 + R / (2 * k)
            #     rMinus = K1 / 4 + K3 * (d_eta * d_eta) / 2 - R / (2 * k)
                
            #     d_force = (wStar + (K1 / (2 * rPlus) * d_eta) + (rMinus / (rPlus)*d_etaPrev)) \
		    #             / (1 / rPlus + (d_J[l]*d_I[l] * kSq) / (rho * A * (1 + sigma_0 * k))) # + (1.0 / hsq * kSq) / (plate_rho * plate_H * (1 + plate_sigma_0 * k)));
                  
            #     force_output[i] += d_force
            #     wNext[d_x - 1 + l] = wStar - (d_J[l]*d_I[l])*((kSq*d_force)/(rho * A * (1+sigma_0*k)))

            d_eta = w[d_x]
            d_etaPrev = wPrev[d_x]
            
            wStar = wNext[d_x]
            
            rPlus = K1 / 4 + K3 * (d_eta * d_eta) / 2 + R / (2 * k)
            rMinus = K1 / 4 + K3 * (d_eta * d_eta) / 2 - R / (2 * k)
                
            d_force = ((wStar + (K1 / (2 * rPlus) * d_eta) + (rMinus / (rPlus)*d_etaPrev)) \
		                / (1 / rPlus + (1/h * kSq) / (rho * A * (1 + sigma_0 * k)))) * 0.5
                  
            force_output[i] += d_force
            wNext[d_x] = wStar - (1/h)*((kSq*d_force)/(rho * A * (1+sigma_0*k)))

        #read output at tip of tine
        output[i] = w[M_w]

        #update state vectors
        uPrev = np.copy(u)
        u = np.copy(uNext)
        
        wPrev = np.copy(w)
        w = np.copy(wNext)
    
    return output, f_out, x_out, force_output
    

output, f_out, x_out, d_force = sim()
plt.plot(x_out, f_out)
plt.show()

e = np.arange(len(output))*1/(48000 * 32)

plt.plot(e, output * 1000)
plt.xlabel('Time (s)')
plt.show()

plt.plot(d_force)
plt.show()