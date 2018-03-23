function transmitTime=...
    findTransTime(sampleNum,readyChnList,svTimeTable,trackResults)
% findTransTime finds the transmitting time of each satellite at a specified
% sample number using the interpolation
% function transmitTime=...
%    findTransTime(sampleNum,readyChnList,svTimeTable,trackResults)
%   Inputs:
%       sampleNum     - absolute sample number from the tracking loop
%       readyChnList  - a list of satellites ready for nav solution
%       svTimeTable   - the transmitting time table
%       trackResults  - output from the tracking function
%
%   Outputs:
%       transmitTime  - transmitting time all ready satellites
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

% Initialize the transmitting time
transmitTime=zeros(1,length(readyChnList));
% size(readyChnList)

% Calculate the transmitting time of each satellite using interpolations
% for channelNr = readyChnList
% %     sampleNum;
%     % Find the index of the sampleNum in the tracking results
%     i = channelNr
%     index_a=max(find((trackResults(channelNr).absoluteSample) <= sampleNum));
%     index_b=min(find((trackResults(channelNr).absoluteSample) >= sampleNum));
%     x1=trackResults(channelNr).absoluteSample(index_a:index_b);
%     y1=[index_a,index_b];
%     index_c=interp1(x1,y1,sampleNum);
%     x2=[index_a,index_b];
%     y2=svTimeTable(channelNr).time(x2);
%     % Find the transmitting time based on the index calculated
%     transmitTime(channelNr)=interp1(x2,y2,index_c);
%     i = channelNr
% end
for channelNr = 1: length(readyChnList)
%     sampleNum;
    % Find the index of the sampleNum in the tracking results
    index_a=max(find((trackResults(channelNr).absoluteSample) <= sampleNum));
    index_b=min(find((trackResults(channelNr).absoluteSample) >= sampleNum));
    x1=trackResults(channelNr).absoluteSample(index_a:index_b);
    y1=[index_a,index_b];
    index_c=interp1(x1,y1,sampleNum);
    x2=[index_a,index_b];
    y2=svTimeTable(channelNr).time(x2);
    % Find the transmitting time based on the index calculated
    transmitTime(channelNr)=interp1(x2,y2,index_c);
end