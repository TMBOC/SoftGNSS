%-----------------------------------------------------------------------------------
% This code has been adapted by Xin Zhang for purposes of course
% "AV423 Satellite Navigation" taught at School of Aeronautics & Astronautics, 
% Shanghai Jiao Tong University,
% from the SoftGNSS v3.0 code base developed for the
% text: "A Software-Defined GPS and Galileo Receiver: A Single-Frequency Approach"
% by Borre, Akos, et.al.
%-----------------------------------------------------------------------------------
function tri = R(tau,T)
tri = (1-abs(tau)./T).*heaviside(1-abs(tau)./T); 
%%%%%%%%%%%%%  R.m  %%%%%%%%%%%%%%%%%%%%

%  HEAVISIDE    Step function.
%      HEAVISIDE(X) is 0 for X < 0, 1 for X > 0, and NaN for X == 0.
%      HEAVISIDE(X) is not a function in the strict sense.
%      See also dirac.

