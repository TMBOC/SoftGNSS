function si = R_BL(tau,b)
% Eq. (2.31) in Winkel (2002)
a = 2*pi*b;
si = (tau+1).*sinint(a.*(tau+1))./pi + cos(a.*(tau+1))./(a*pi) + (tau-1).*sinint(a.*(tau-1))./pi ...
            +cos(a.*(tau-1))./(a*pi) - 2.*tau.*sinint(a.*tau)./pi - cos(a.*tau)./(pi*pi*b); 
%%%%%%%%%%%%%  R_BL.m  %%%%%%%%%%%%%%%%%%%%
