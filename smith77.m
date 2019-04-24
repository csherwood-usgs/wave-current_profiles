% smith77 - Analytical solution for wbl with linear eddy viscosity
% Smith, J.Dungan (1977) Chap. 13 in The Sea, vol. 6
u0 = 100;
% w = 2*pi/8;
w = 1.59;

% Ko = vk*ustrm;
Ko = 2.2;
zo = 1e-3;
z = logspace(log10(zo),log10(100),50)';
nz = length(z);

x = 2*sqrt(w*z/Ko);
x0 = 2*sqrt(w*zo/Ko);

% See Abromwitz & Stegun, 9.9.2, p. 379
bex = besselk(0,x.*exp(pi*1i/4));
kerx = real(bex);
keix = imag(bex);
bex0 = besselk(0,x0.*exp(pi*1i/4));
kerx0 = real(bex0);
keix0 = imag(bex0);

wt = (0:pi/10:pi);
nt = length(wt);
%%
for i=1:nt
u2(:,i) =u0*( cos(wt(i))...
   -( (kerx*kerx0+keix*keix0)/(kerx0.^2+keix0.^2) ).*cos(wt(i))...
   -( (kerx*keix0-keix*kerx0)/(kerx0.^2+keix0.^2) ).*sin(wt(i)));
end
figure(1);clf
semilogy(u2,z)
title('Smith, 1977, Fig. 3')