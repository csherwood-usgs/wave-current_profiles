% M94_orbital - Calculate velocity profiles using Madsen (1994) w-c bbl model
%
% This model is based on two unconnected eddy viscosity profiles, one inside the wave bbl,
% and one above it...so it produces a discontinuous velocity profile.
%
% Calls:
%   m94r - Solves the model
%   ker - real part of first-order Kelvin function
%   kei - imag part of first-order Kelvin function
%
% Madsen, O.S. (1994). Spectral Wave-Current Bottom Boundary Layer Flows.
% Coastal Engineering 1994. https://doi.org/10.1061/9780784400890.030
%
% csherwood@usgs.gov

% Enter values here:
ucr = .2 % current speed at zr [m/s]
ubr = .5 % representative wave-orbital velocity [m/s]
T = 10   % wave period [s]
phiwc=0. % angle between waves and currents [radians]
zr = 1.  % reference elevation where currents are observed [m]
kN = 1e-3 % bottom Nikuradse roughness (~ grain size) [m]

% Begin calcs
vk = 0.41 % von Karman's constant
zo = kN/30
wr = 2*pi/T % angular frequency [s-1]
iverbose = 0;

% Madsen 94 w-c bbl model
m = m94r( ubr, wr, ucr, zr, phiwc, kN, iverbose )

% elevations for profiles
nz = 20;
z = logspace(log10(zo), log10(zr), nz)'; % [m]
delw = m.dwc        % M94 eqn. 36 wave boundary layer thickness [m]
zetan = z./delw;    % M94 eqn. 8
zetano = zetan(1);

% orbital velocity time series (radians)
nt = 10;
wt = ( 0: pi/(nt-1): pi);

% wave-orbital velocity profiles
ubp = zeros(nz,nt);
for j=1:nt
   num = ( ker(2*sqrt(zetan))  + 1i*kei(2*sqrt(zetan)))./ ...
         ( ker(2*sqrt(zetano)) + 1i*kei(2*sqrt(zetano)));
   ubp(:,j) = ubr*real( (1.-num) * real(exp(1i*wt(j))) ); % M94 eqn. 6
end

% mean current profile [m/s] - x is in wave direction, y is perpendicular
% M94 eqns. 9 and 10
ucx = zeros(nz,1);
ucy = zeros(nz,1);
below = (z<=delw);
ucx(z<=delw) = (m.ustrc.^2/(vk*m.ustrr))*log(z(z<=delw)./zo) .*cos(phiwc);
ucx(z>delw)  = (m.ustrc/vk)*log(z(z>delw)./m.zoa) .*cos(phiwc);
ucy(z<=delw) = (m.ustrc.^2/(vk*m.ustrr))*log(z(z<=delw)./zo) .*sin(phiwc);
ucy(z>delw)  = (m.ustrc/vk)*log(z(z>delw)./m.zoa) .*sin(phiwc);

% plots
figure(1); clf
plot(1.2*[-ubr ubr],[m.dwc,m.dwc],'--b')
hold on
plot(ubp,z)
plot(ucx,z,'-b','linewidth',2)
plot(ucy,z,'-r','linewidth',2)
set(gca,'yscale','log')
ylabel('Elevation [m]')
xlabel('Velocity [m/s]')
title('Wave Orbital and Mean Current Profiles')
print -dpng profs.png

figure(2); clf
plot(1.2*[-ubr ubr],[m.dwc,m.dwc],'--b')
hold on
h1=plot(ubp(:,1)+ucx,z,'-b')
hold on
h2=plot(ucy,z,'-r')
set(gca,'yscale','log')
ylabel('Elevation [m]')
xlabel('Velocity [m/s]')
legend([h1;h2],'parallel to wave direction','perpendicular to wave direction')
title('Max. Combined Velocity Profiles')
print -dpng comb_profs.png