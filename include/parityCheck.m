function word = parityCheck(word, D29Star, D30Star)
%Checks the parity of the supplied 30bit word.If the parity check is not
%passed, the program will be quit with error messages shown
%The last two parity bits of the previous word is used for the calculation.
%A note on the procedure is supplied by the GPS standard positioning
%service signal specification.
%
%word = parityCheck(word, D29Star, D30Star)
%
%   Inputs:
%       word        - an array with 30 bit long word from the navigation
%                   message (a character array, must contain only '0' or
%                   '1').
%       D29Star     - the 29th bit of the previous word (char type).
%       D30Star     - the 30th bit of the previous word (char type).
%
%   Outputs:
%       word        - word with corrected polarity of the data bits
%                   (character array).

%--------------------------------------------------------------------------
% Copyright (C) D.M.Akos
% Written by Xiaofan Li
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

% Convert the word,D29Star,D30Start from 0,1 
% in char type to 1,-1 in double type 
word=1-2*(word-48);
D29Star=1-2*(D29Star-48);
D30Star=1-2*(D30Star-48);

% Correct the polarity of the data word if necessary
word(1:24)=word(1:24)*D30Star;

% Establish matrices for the parity check of the last six bits
h1=[1 1 1 0 1 1 0 0 0 1 1 1 1 1 0 0 1 1 0 1 0 0 1 0];
h2=[0 1 1 1 0 1 1 0 0 0 1 1 1 1 1 0 0 1 1 0 1 0 0 1];
h3=[1 0 1 1 1 0 1 1 0 0 0 1 1 1 1 1 0 0 1 1 0 1 0 0];
h4=[0 1 0 1 1 1 0 1 1 0 0 0 1 1 1 1 1 0 0 1 1 0 1 0];
h5=[1 0 1 0 1 1 1 0 1 1 0 0 0 1 1 1 1 1 0 0 1 1 0 1];
h6=[0 0 1 0 1 1 0 1 1 1 1 0 1 0 1 0 0 0 1 0 0 1 1 1];
H=[h1;h2;h3;h4;h5;h6];
d=word(1:24);
DStarMat=[D29Star,D30Star,D29Star,D30Star,D30Star,D29Star];

% Calculation of the D
for i=1:6
    temp=H(i,:).*d;
    D(i)=prod(temp(find(temp~=0)))*DStarMat(i);
end;

% Convert the word and D back into 0,1 in char type
word=dec2bin((1-word)/2);
D=dec2bin((1-D)/2);

% If the D does not match the last six bits in the word
% An error message will be reported
if D~=word(25:30)
    error('Parity Check Failed!');
end

