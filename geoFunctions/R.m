function tri = R(tau,T)
tri = (1-abs(tau)./T).*heaviside(1-abs(tau)./T); 
%%%%%%%%%%%%%  R.m  %%%%%%%%%%%%%%%%%%%%

%  HEAVISIDE    Step function.
%      HEAVISIDE(X) is 0 for X < 0, 1 for X > 0, and NaN for X == 0.
%      HEAVISIDE(X) is not a function in the strict sense.
%      See also dirac.

