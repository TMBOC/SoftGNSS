function transmitTime=transTimeTable(activeChnList,trackResults...
    ,subFrameStart,TOW,settings)
%transTimeTable establishes a transmitting time table for all satellites on
%track. TOW is the transmiting time at the start of the subframe for all
%satellites, based on that the time table can be established.
%transmitTime=transTimeTable(activeChnList,trackResults...
%    ,subFrameStart,TOW,settings)
%
%   Inputs:
%       activeChnList - a list of active satellites
%       trackResults  - output from the tracking function
%       subFrameStart - the start of first detected subframe
%       TOW           - transmitting time of satellites at the first
%                       detected subframe
%       settings      - receiver settings
%
%   Outputs:
%       transmitTime  - Transmitting time table for all active satellites
%--------------------------------------------------------------------------
% Copyright (C) D.M.Akos
% Written by Xiaofan Li at University of Colorado at Boulder
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

% Initialize the transmitting time table
transmitTime.time=zeros(1,settings.msToProcess);
transmitTime.PRN=[];
transmitTime = repmat(transmitTime, 1, max(activeChnList));

% Establish the time table based on the TOW and its position
for channelNr = activeChnList
    transmitTime(channelNr).PRN=trackResults(channelNr).PRN;
    for i=1:settings.msToProcess
        transmitTime(channelNr).time(i)=...
            TOW-subFrameStart(channelNr)*0.001+i*0.001;
    end;
end