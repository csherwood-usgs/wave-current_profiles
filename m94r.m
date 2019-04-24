function m = m94r( ubr, wr, ucr, zr, phiwc, kN, iverbose )
% M94 - Grant-Madsen model from Madsen(1994)
% function m =  m94( ubr, wr, ucr, zr, phiwc, kN, iverbose )
%
%   Input:
%      ubr = rep. wave-orbital velocity amplitude outside wbl [m/s]
%      wr = rep. angular wave frequency = 2pi/T [rad/s]
%      ucr = current velocity at height zr [m/s]
%      zr = reference height for current velocity [m]
%      phiwc = angle between currents and waves at zr (radians)
%      kN = bottom roughness height (e.q. Nikuradse k) [m]
%      iverbose = switch; when 1, extra output
%   Returned in structure m:
%      m.ustrc  = current friction velocity         u*c [m/s]
%      m.ustrr  = w-c combined friction velocity    u*r [m/s]
%      m.ustrwm = wave max. friction velocity      u*wm [m/s]
%      m.dwc = wave boundary layer thickness [m]
%      m.fwc = wave friction factor [ ]
%      m.zoa = apparent bottom roughness [m]
%

% Chris Sherwood, USGS
% Last revised December 17, 2001

% This version prints diagnostics to help with constraining parameter space

MAXIT = 20;
PI = 3.1415926535;
VK = 0.41;
rmu=zeros(MAXIT,1);
Cmu=zeros(MAXIT,1);
fwci=zeros(MAXIT,1);
dwci=zeros(MAXIT,1);
ustrwm2=zeros(MAXIT,1);
ustrr2=zeros(MAXIT,1);
ustrci=zeros(MAXIT,1);

%     ...junk return values
ustrc = 99.99;
ustrwm = 99.99;
ustrr = 99.99;
fwc = .4;
kN = min(kN, 0.9*zr); % overall limit on roughness
zoa = kN/30.;
dwc = kN; 
phicwc = phiwc;

% ...some data sanity checks
if( wr <= 0. )
   fprintf(1,'WARNING: Bad value for frequency in M94: wr=%g\n',...
      wr );
   return;
end
if( ubr < 0. )
   fprintf(1,'WARNING: Bad value for orbital vel. in M94: ub=%g\n',...
      ubr);
   return;
end
if( kN < 0. )
   fprintf(1,'WARNING: Wierd value for roughness in M94: kN=%g'\n',...
      kN);
   return;
end
if( (zr<zoa || zr<0.05) && iverbose==1)
   fprintf(1,'WARNING: Low value for ref. level in M94: zr=%g\n',...
      zr);
end

zo = kN/30.;
if(ubr <= 0.01),
   if(ucr <= 0.01),
      %     ...no waves or currents
      ustrc = 0.;
      ustrwm = 0.;
      ustrr = 0.;
      m.ustrc = ustrc;
      m.ustrr = ustrr;
      m.ustrwm = ustrwm;
      m.dwc = dwc;
      m.fwc = fwc;
      m.zoa = zoa;
      return
   end
   %  ...no waves
   ustrc = ucr * VK / log(zr/zo);
   ustrwm = 0.;
   ustrr = ustrc;
   m.ustrc = ustrc;
   m.ustrr = ustrr;
   m.ustrwm = ustrwm;
   m.dwc = dwc;
   m.fwc = fwc;
   m.zoa = zoa;
   return
end

cosphiwc =  abs(cos(phiwc));
rmu(1) = 0.;
Cmu(1) = 1.;
fwci(1) = fwc94( Cmu(1), (Cmu(1)*ubr/(kN*wr)) ); %Eqn. 32 or 33
ustrwm2(1)= 0.5*fwci(1)*ubr*ubr;                 %Eqn. 29
ustrr2(1) = Cmu(1)*ustrwm2(1);                   %Eqn. 26
ustrr = sqrt( ustrr2(1) );
dwci(1) = kN;
if ((Cmu(1)*ubr/(kN*wr)) >= 8.), dwci(1)= min( 2.*VK*ustrr/wr, 0.9*zr); end
lnzr = log(zr/dwci(1));
lndw = log(dwci(1)/zo);
lnln = lnzr/lndw;
bigsqr = (-1.+sqrt(1+ ((4.*VK*lndw)/(lnzr*lnzr))*ucr/ustrr ));
ustrci(1) = 0.5*ustrr*lnln*bigsqr;
nit = 1;
diffw = 1.;
for i=2:MAXIT
   rmu(i) = ustrci(i-1)*ustrci(i-1)/ustrwm2(i-1);
   Cmu(i) = sqrt(1.+2.*rmu(i)*cosphiwc+rmu(i)*rmu(i));%Eqn 27
   fwci(i) = fwc94( Cmu(i), (Cmu(i)*ubr/(kN*wr)) );   %Eqn. 32 or 33
   ustrwm2(i)= 0.5*fwci(i)*ubr*ubr;                   %Eqn. 29
   ustrr2(i) = Cmu(i)*ustrwm2(i);                     %Eqn. 26
   ustrr = sqrt( ustrr2(i) );
   dwci(i) = kN;
   if ((Cmu(i)*ubr/(kN*wr))>= 8.) dwci(i)= min( 2.*VK*ustrr/wr, 0.9*zr); end % Eqn.36...also needs to be limited
   lnzr = log( zr/dwci(i) );
   lndw = log( dwci(i)/zo );
   lnln = lnzr/lndw;
   bigsqr = (-1.+sqrt(1+ ((4.*VK*lndw)/(lnzr*lnzr))*ucr/ustrr ));
   ustrci(i) = 0.5*ustrr*lnln*bigsqr;                  %Eqn. 38
   diffw = abs( (fwci(i)-fwci(i-1))/fwci(i) );
   if(diffw < 0.000005), break, end
   nit = nit+1;
end
ustrwm = sqrt( ustrwm2(nit) );
ustrc = ustrci(nit);
ustrr = sqrt( ustrr2(nit) );
phicwc = phiwc;
zoa = exp( log(dwci(nit))-(ustrc/ustrr)*log(dwci(nit)/zo) ); %Eqn. 11
fwc = fwci(nit);
dwc = dwci(nit);
if(iverbose==1)
   for i=1:nit
      fprintf(1,...
         'i=%2d fwc=%9.6f kN=%9.6f dwc=%9.6f zoa=%9.6f u*c=%9.4f u*wm=%9.4f u*r=%9.4f\n',...
         i,fwci(i),kN, dwci(i),zoa, ustrci(i),sqrt(ustrwm2(i)),...
         sqrt(ustrr2(i)))
   end
end
if(iverbose==2)
   i=nit;
   fprintf(1,...
      'i=%2d fwc=%9.6f kN=%9.6f dwc=%9.6f zoa=%9.6f u*c=%9.4f u*wm=%9.4f u*r=%9.4f\n',...
      i,fwci(i),kN, dwci(i), zoa, ustrci(i), sqrt(ustrwm2(i)),...
      sqrt(ustrr2(i)))
end
m.ustrc = ustrc;
m.ustrr = ustrr;
m.ustrwm = ustrwm;
m.dwc = dwc;
m.fwc = fwc;
m.zoa = zoa;
return

function fwc = fwc94( cmu, cukw )
% FWC94r - Wave-current friction factor
% Equations 32 and 33 in Madsen, 1994

% alternative to next two lines
% fwc = .43; % value at cukw = 0.2 ...this is pretty rough!
% if( (cukw >= 0.2) && (cukw <= 100.) )

fwc = .3; % value at cukw = 0.352...a little less rough
if( (cukw >= 0.352) && (cukw <= 100.) )
   fwc = cmu*exp( 7.02*cukw^(-0.078)-8.82 );
elseif( cukw>100. )
   fwc = cmu*exp( 5.61*cukw^(-0.109)-7.30 );
end
