function ki = kei(x)
% kei - imag part of Kelvin function of zero order
% Abromowitz & Stegun, 9.9.2, p. 379
ki = imag(besselk(0,x.*exp(pi*1i/4)));