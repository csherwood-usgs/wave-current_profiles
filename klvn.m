function k = klvn(x)
% Zero order Kelvin functions of x
% Polynomial approximation from Abromowitz and Stegun,
% Handbook of Mathematical Functions, p. 384.
%
% ber and bei valid for -8 < x 8
% ker and kei valid for 0 < x < 8
%
% Adpated from KLVN.F in NBM87 by Scott Glenn
%
% csherwood@usgs.gov Nov 2009
% Note that in Matlab (see p. 379)
% ber + 1i*bei = besselj(0,x.*exp(3*pi*1i/4))
% ker + 1i*kei = besselk(0,x.*exp(pi*1i/4))

x2=x./2;
x8=x./8;

ber=1-64*(x8.^4)+113.77777774*(x8.^8)-32.36345652*(x8.^12)...
   +2.64191397*(x8.^16)-.08349609*(x8.^20)+.00122552*(x8.^24)...
   -.00000901*(x8.^28);

bei=16*(x8.^2)-113.77777774*(x8.^6)+72.81777742*(x8.^10)...
   -10.56765779*(x8.^14)+.52185615*(x8.^18)-.01103667*(x8.^22)...
   +.00011346*(x8.^26);

ker=-ber.*log(x2)+pi/4 .*bei-.57721566-59.05819744*(x8.^4)...
   +171.36272133*(x8.^8)-60.60977451*(x8.^12)...
   +5.65539121*(x8.^16)-.19636347*(x8.^20)...
   +.00309699*(x8.^24)-.00002458*(x8.^28);

kei=-bei.*log(x2)-pi/4 .*ber+6.76454936*(x8.^2)...
   -142.91827687*(x8.^6)+124.23569650*(x8.^10)...
   -21.30060904*(x8.^14)+1.17509064*(x8.^18)...
   -.02695875*(x8.^22)+.00029532*(x8.^26);

% 9.10.31 on p. 383
xinv = 1 ./x;
sq2 = sqrt(2);
ker2kei2 = ( 0.5*pi*xinv ).*exp(-x*sq2) .* ...
   (1 - xinv./(4*sq2).*xinv + (xinv.^2)./64 + ...
   (xinv.^3).*(33/(256*sq2))-(xinv.^4).*(1797/8192));
k.x = x;
k.ber = ber;
k.bei = bei;
k.ker = ker;
k.kei = kei;
k.ker2kei2 = ker2kei2;

return