function [CNo]= CNoVSM(I,Q,T)
%function [CNo]= CNoVSM(I,Q,T)
%Calculate CNo using the Variance Summing Method
%
%[CNo]= CNoVSM(I,Q)
%
%   Inputs:
%       I           - Prompt In Phase values of the signal from Tracking
%       Q           - Prompt Quadrature Phase values of the signal from Tracking
%       T          - Accumulation interval in Tracking (in sec)
%   Outputs:
%       CNo         - Estimated C/No for the given values of I and Q
%
%
%
%--------------------------------------------------------------------------
% Copyright (C) D.M.Akos
% Written by Sirish Jetti
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------


%Calculate Power
Z=I.^2+Q.^2;
%Calculate the mean and variance of the Power
Zm=mean(Z);
Zv=var(Z);
%Calculate the average carrier power
Pav=sqrt(Zm^2-Zv);
%Calculate the variance of the noise
Nv=0.5*(Zm-Pav);
%Calculate C/No
CNo=10*log10(abs((1/T)*Pav/(2*Nv)));