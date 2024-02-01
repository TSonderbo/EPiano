import scipy as sp
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from numba import njit
    
fs = 48000 #sample rate
k = 1/fs #sampling Period

lenSound = fs #length of sound in samples

L = 1
rho = 7850 #material Density
r = 5 * 10**-4
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
uNext = np.zeros((N+1,1))
u = np.zeros((N+1,1))
vg1Prev = 0

# raised cosine
loc = round(0.8 * N); # Center location
halfWidth = round(N/10); # Half-width of raised cosine
width = 2 * halfWidth; # Full width
rcX = range(0,width); # x-locations for raised cosine
rc = np.array([0.5 - 0.5 * math.cos(2 * math.pi * x / width) for x in rcX]) # raised cosine
u[loc-halfWidth : loc+halfWidth] = rc.reshape(-1,1) # initialise current state

#output
outLoc = N #ouput location
output = np.zeros((lenSound,1)) #output vector

uPrev = np.copy(u)

vg1Prev = 0
    
def animate(i):
    
    global u
    global uPrev
    global uNext
    global vg1Prev
    
    for j in range (0,6):
        for l in range(2, N-1):
            uNext[l] = ((2-6*muSq-((4*sigma_1*k)/hSq))*u[l] \
                + (4*muSq + ((2*sigma_1*k)/hSq))*(u[l+1] + u[l-1]) \
                - muSq*(u[l+2] + u[l-2]) + (-1+sigma_0*k+((4*sigma_1*k)/hSq))*uPrev[l] \
                - ((2*sigma_1*k)/hSq) * (uPrev[l+1] + uPrev[l-1])) \
                /(1+sigma_0*k)   
        
        vg1 = 2*u[N] - u[N-1]
        vg2 = u[N-2]-2*u[N-1]+2*vg1

        l = N-1
        uNext[l] = ((2-6*muSq-(4*sigma_1*k)/hSq)*u[l] \
            + (4*muSq + (2*sigma_1*k)/hSq)*(u[l+1] + u[l-1]) \
            - muSq*(vg1 + u[l-2]) + (-1+sigma_0*k+(4*sigma_1*k)/hSq)*uPrev[l] \
            - ((2*sigma_1*k)/hSq) * (uPrev[l+1] + uPrev[l-1]))/(1+sigma_0*k)   

        l = N
        uNext[l] = ((2-6*muSq-(4*sigma_1*k)/hSq)*u[l] \
            + (4*muSq + (2*sigma_1*k)/hSq)*(vg1 + u[l-1]) \
            - muSq*(vg2 + u[l-2]) + (-1+sigma_0*k+(4*sigma_1*k)/hSq)*uPrev[l] \
            - ((2*sigma_1*k)/hSq) * (vg1Prev + uPrev[l-1]))/(1+sigma_0*k)   

        vg1Prev = vg1
        
        uPrev = np.copy(u)
        u = np.copy(uNext)
    
    line.set_ydata(u)
    
    return line,
    
fig, ax = plt.subplots()
    
x = np.arange(N+1)

line = ax.plot(x,u)[0]

ax.set_ylim((-1, 1))
    
# calling the animation function      
anim = animation.FuncAnimation(fig, animate, frames = 6000, interval = 1, blit = True)  
   
# saves the animation in our desktop 
anim.save('bar.mp4', writer = 'ffmpeg', fps = 30) 