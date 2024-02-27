clear all;

%% Initialise variables
fs = 48000*4; % Sample rate [Hz]
k = 1 / fs; % Time step [s]

lengthSound = fs; % Length of the simulation (1 second) [samples]

L = 0.1564; % Length [m]
r = 2 * 10 ^-3; % Radius
A = pi * r^2; % Cross-sectional Area
rho = 7850; %Material Density
E = 2 * 10^11; %Young's Modulus
I = pi * r^4 / 4; %Inertia

sigma_0 = 1; %Frequency independent damping
sigma_1 = 0.05; %Frequency dependent damping

kappa = sqrt((E * I) / (rho * A)); %Stiffness coefficient

h = sqrt((4 * sigma_1 * k + sqrt((4 * sigma_1 * k)^2 + 16 * kappa * kappa * k * k)) * 0.5); %Grid spacing [m]
N = floor(L/h); %Number of intervals between grid points
h = L / N; %Recalculation of grid spacing based on integer N

mu = (kappa * k)/ h^2;
muSq = mu^2;
hSq = h^2;

%Initialise state vectors
uNext = zeros(N+1, 1);
u = zeros(N+1, 1);

%Raised cosine
loc = round(0.8 * N); %Center location
halfWidth = round(N/10); %Half-width of raised cosine
width = 2 * halfWidth; %Full width
rcX = 0:width; %x-locations for raised cosine
rc = 0.5 - 0.5 * cos(2 * pi * rcX / width); %Raised cosine
u(loc-halfWidth : loc+halfWidth) = rc; %Initialise current state

%Set initial velocity to zero
uPrev = u;

%Range of calculation
range = 1+2:N-1;

%Virtual grid points
vg1 = 0;
vg2 = 0;
vg1Prev = 0;

%==================================================================
%Simulation loop
%==================================================================

outLoc = N+1;
out = zeros(lengthSound, 1);

for n = 1:lengthSound

    %Update equation
    uNext(range) = ((2 - 6*muSq - ((4 * sigma_1 * k) / hSq)) * u(range) ...
			+ (4 * muSq + ((2 * sigma_1 * k) / hSq)) * (u(range + 1) + u(range - 1)) ...
			- muSq * (u(range + 2) + u(range - 2)) + (-1 + sigma_0 * k + ((4 * sigma_1 * k) / hSq)) * uPrev(range) ...
			- ((2 * sigma_1 * k) / hSq) * (uPrev(range + 1) + uPrev(range - 1))) ...
			/(1+sigma_0*k);

    %Update virtual grid points
    vg1 = 2 * u(N+1) - u(N);
    vg2 = u(N-1) - 2*u(N) + 2*vg1;

    %Update Free Boundary
    uNext(N) = ((2 - 6*muSq - ((4 * sigma_1 * k) / hSq)) * u(N) ...
			+ (4 * muSq + ((2 * sigma_1 * k) / hSq)) * (u(N + 1) + u(N - 1)) ...
			- muSq * (vg1 + u(N - 2)) + (-1 + sigma_0 * k + ((4 * sigma_1 * k) / hSq)) * uPrev(N) ...
			- ((2 * sigma_1 * k) / hSq) * (uPrev(N + 1) + uPrev(N - 1))) ...
			/(1+sigma_0*k);

    uNext(N+1) = ((2 - 6*muSq - ((4 * sigma_1 * k) / hSq)) * u(N + 1) ...
			+ (4 * muSq + ((2 * sigma_1 * k) / hSq)) * (vg1 + u(N)) ...
			- muSq * (vg2 + u(N - 1)) + (-1 + sigma_0 * k + ((4 * sigma_1 * k) / hSq)) * uPrev(N+1) ...
			- ((2 * sigma_1 * k) / hSq) * (vg1Prev + uPrev(N))) ...
			/(1+sigma_0*k);

    %Read output
    out(n) = u(outLoc);
    
    %Update states
    vg1Prev = vg1;
    uPrev = u;
    u = uNext;
end

plot(out);
