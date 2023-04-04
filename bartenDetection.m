function output = bartenDetection(u,e,L,D,k,eta0,sigma0,eg,u00)

% Modified Barten CSF model for detection of luminance patterns, proposed by Bozorgian et al. (2022)
% u: Frequency range in cpd
% e: Eccentricity in degrees
% L: Average luminance of the observed object in cd/m^2
% D: Field diameter in degrees
% k: Signal to Noise ratio
% eta0: Constant for quantom efficiency
% sigma0: Constant for the eye MTF
% eg: Eccentricity constant (can be different for various subjects)
% u00: Spatial frequency above which the lateral inhibition ceases as a function of eccentricity in Fovea

% Density of total Retinal ganglion cells in fovea
Nc0       = 12000;

% Density of parasol Retinal ganglion cells in fovea
Ngp0      = 0.05 * 36000;

% Density of parasol Retinal ganglion cells as a function of eccentricity
Ngp       = Ngp0 * (0.85/(1+(e/0.45)^2) + 0.15/(1+(e/eg)^2));

% constants for calculation of sigma based on pupil diameter
Cab       = 0.08;

% Eye Integration Time
T         = 0.1;

% maximum angular size of the integration area of the noise as a function of eccentricity
Xmax0     = 12;
Xmax      = Xmax0 * ( 0.85/(1+(e/4)^2) + 0.15/(1+(e/12)^2) )^-0.5;
Ymax      = Xmax;

% maximum number of cycles over which the eye can integrate the information
Nmax      = 15;

% quantum efficiency of the eye as a function of eccentricity
eta       = eta0 * (0.4/(1+(e/7)^2) + 0.48/(1+(e/20)^2) + 0.12);

% Spectral density of the neural noise as a function of eccentricity
Phi00     = 3*10^-8;
Phi0      = Phi00 * (1+0.24*e);

% Spatial frequency above which the lateral inhibition ceases as a function of eccentricity
u0        = u00 * (Ngp/Ngp0)^0.5 * ( 0.85/(1+(e/4)^2) + 0.13/(1+(e/20)^2) + 0.02 )^-0.5;

% Photon conversion factor
p         = 1.240 *10^6;

% Angular field size of object in degrees
X0        = sqrt(pi/4 * D^2);
Y0        = X0;

% Pupil diameter
d         = 5 - 3 * tanh(0.4 * log10(L*(X0^2/(40^2))));

% standard deviation of the line-spread function caused by the discrete sturcture of the retina (3600 converts degree to arcmin)
rc       = (0.45/30)*e + 0.25;
sigmaret = rc/sqrt(18/5);

% standard deviation of the line-spread function caused by the other parts of eye
sigma00   = sqrt(sigma0^2 - 0.1318^2);

% Standard deviation of the line-spread function resulting from the convolution  of the different elements of convolution process
sigma     = sqrt(sigmaret^2 + sigma00^2 + (Cab*d)^2);

% Moduation transfer function of the optics of the eye 
Mopt      = exp(-2 * pi^2 *sigma^2 .* u.^2 .* (1/3600));

% Modulation transfer function of the lateral inhibition process
Mlat      = sqrt(1 - exp(-(u./u0).^2));

% Retinal illuminance in Troland
E         = (pi*(d^2)/4) .* L .* ( 1-(d/9.7)^2 + (d/12.4)^4 );

% Spectral Density of photon noise
PhiPhoton = 1/(eta*p*E);

% Spatial equations
X         = ( 1/(X0^2) + 1/(Xmax^2) + ( ((0.5*X0)^2 + 4*e^2)/((0.5*X0)^2) + e^2) * (u.^2)./(Nmax^2) ).^-0.5;
Y         = ( 1/(Y0^2) + 1/(Ymax^2) + ( ((0.5*X0)^2)/((0.5*X0)^2) + e^2) * (u.^2)./(Nmax^2) ).^-0.5;

% Contrast Sensitivity Function
% add * sqrt(2) in case of binocular mode 
output    = (Mopt./(2*k)) .* sqrt( (X.*Y.*T) ./ (PhiPhoton + Phi0./(Mlat.^2)));

end
