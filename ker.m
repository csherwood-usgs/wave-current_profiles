function kr = ker(x)
% ker - real part of Kelvin function of zero order
% Abromowitz & Stegun, 9.9.2, p. 379
kr = real(besselk(0,x.*exp(pi*1i/4)));