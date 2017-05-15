function [navSolutions, eph,svTimeTable,activeChnList] = postNavigation(trackResults, settings)
%Function calculates navigation solutions for the receiver (pseudoranges,
%positions). At the end it converts coordinates from the WGS84 system to
%the UTM, geocentric or any additional coordinate system.
%
%[navSolutions, eph] = postNavigation(trackResults, settings)
%
%   Inputs:
%       trackResults    - results from the tracking function (structure
%                       array).
%       settings        - receiver settings.
%   Outputs:
%       navSolutions    - contains measured pseudoranges, receiver
%                       clock error, receiver coordinates in several
%                       coordinate systems (at least ECEF and UTM).
%       eph             - received ephemerides of all SV (structure array).

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%
% Copyright (C) Darius Plausinaitis
% Written by Darius Plausinaitis with help from Kristin Larson
% Modified By Xiaofan Li at University of Colorado at Boulder
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

%CVS record:
%$Id: postNavigation.m,v 1.1.2.22 2006/08/09 17:20:11 dpl Exp $


%% Check is there enough data to obtain any navigation solution ===========
% It is necessary to have at least three subframes (number 1, 2 and 3) to
% find satellite coordinates. Then receiver position can be found too.
% The function requires all 5 subframes, because the tracking starts at
% arbitrary point. Therefore the first received subframes can be any three
% from the 5.
% One subframe length is 6 seconds, therefore we need at least 30 sec long
% record (5 * 6 = 30 sec = 30000ms). We add extra seconds for the cases,
% when tracking has started in a middle of a subframe.
tic
if (settings.msToProcess < 36000) || (sum([trackResults.status] ~= '-') < 4)
    % Show the error message and exit
    disp('Record is to short or too few satellites tracked. Exiting!');
    navSolutions = [];
    eph          = [];
    svTimeTable  = [];
    activeChnList = [];
    return
end

%% Find preamble start positions ==========================================

[subFrameStart, activeChnList] = findPreambles(trackResults, settings);
%% Decode ephemerides =====================================================

for channelNr = activeChnList
    
    %=== Convert tracking output to navigation bits =======================
    
    %--- Copy 5 sub-frames long record from tracking output ---------------
    navBitsSamples = trackResults(channelNr).I_P(subFrameStart(channelNr) - 40 : ...
        subFrameStart(channelNr) + (1500 * 20) -1)';
    
    %--- Group every 20 vales of bits into columns ------------------------
    navBitsSamples = reshape(navBitsSamples, ...
        20, (size(navBitsSamples, 1) / 20));
    
    %--- Sum all samples in the bits to get the best estimate -------------
    navBits = sum(navBitsSamples);
    
    %--- Now threshold and make 1 and 0 -----------------------------------
    % The expression (navBits > 0) returns an array with elements set to 1
    % if the condition is met and set to 0 if it is not met.
    navBits = (navBits > 0);
    
    %--- Convert from decimal to binary -----------------------------------
    % The function ephemeris expects input in binary form. In Matlab it is
    % a string array containing only "0" and "1" characters.
    navBitsBin = dec2bin(navBits);
    %=== Decode ephemerides and TOW of the first sub-frame ================
    [eph(trackResults(channelNr).PRN), TOW] = ...
        ephemeris(navBitsBin(3:1502)', navBitsBin(1),navBitsBin(2));
%     old version of ephemeris.m
%     [eph(trackResults(channelNr).PRN), TOW] = ...
%                             ephemeris(navBitsBin(3:1502)', navBitsBin(1));
    
    %--- Exclude satellite if it does not have the necessary nav data -----
    % If the satellite accuracy or health is not in reliable values, then
    % this satellite is excluded as well
    if (isempty(eph(trackResults(channelNr).PRN).IODC) || ...
            isempty(eph(trackResults(channelNr).PRN).IODE_sf2) || ...
            isempty(eph(trackResults(channelNr).PRN).IODE_sf3) || ...
            eph(trackResults(channelNr).PRN).accuracy >=3 ||...
            eph(trackResults(channelNr).PRN).health~=0)
        
        %--- Exclude channel from the list (from further processing) ------
        activeChnList = setdiff(activeChnList, channelNr);
        s=sprintf('PRN %d is excluded from the active channel',...
            trackResults(channelNr).PRN);
        disp(s);
    end
end

%% Check if the number of satellites is still above 3 =====================
if (isempty(activeChnList) || (size(activeChnList, 2) < 4))
    % Show error message and exit
    disp('Too few satellites with ephemeris data for postion calculations. Exiting!');
    navSolutions = [];
    eph          = [];
    svTimeTable  = [];
    activeChnList = [];
    return
end

%% Initialization =========================================================
% Set the satellite elevations array to INF to include all satellites for
% the first calculation of receiver position. There is no reference point
% to find the elevation angle as there is no receiver position estimate at
% this point.
satElev  = inf(1, settings.numberOfChannels);

% Save the active channel list. The list contains satellites that are
% tracked and have the required ephemeris data. In the next step the list
% will depend on each satellite's elevation angle, which will change over
% time.
readyChnList = activeChnList;
% Establish the transmitting time table
svTimeTable.time=zeros(1,settings.msToProcess);
svTimeTable.PRN=[];
svTimeTable = repmat(svTimeTable, 1, max(activeChnList));

% Establish the time table based on the TOW and its position
for channelNr = activeChnList
    svTimeTable(channelNr).PRN=trackResults(channelNr).PRN;
    for i=1:settings.msToProcess
        svTimeTable(channelNr).time(i)=...
            TOW-subFrameStart(channelNr)*0.001+(i-1)*0.001;
    end;
end
% 
% svTimeTable=...
%     transTimeTable(activeChnList,trackResults,subFrameStart,TOW,settings);
% transmitTime = TOW;

% Find the last sample number in the tracking results
lastSample=inf(1,max(readyChnList));

for channelNr = readyChnList
    lastSample(channelNr) = ...
        trackResults(channelNr).absoluteSample(end);
end
% Find the step size for navigation solution
navStep=settings.samplingFreq/settings.navSolRate;

%##########################################################################
%#   Do the satellite and receiver position calculations                  #
%##########################################################################

%% Initialization of current measurement ==================================
for currMeasNr =1:fix((min(lastSample) - ...
        settings.samplingFreq/settings.navSolRate-...
        settings.skipNumberOfSamples) /navStep);
    currMeasNr
% for currMeasNr =1:fix((settings.msToProcess - max(subFrameStart)) / ...
%                                                      (1000/settings.navSolRate))
    % Exclude satellites, that are below elevation mask
    activeChnList = intersect(find(satElev >= settings.elevationMask), ...
        readyChnList);
    
    % Save list of satellites used for position calculation
    navSolutions.channel.PRN(activeChnList, currMeasNr) = ...
        [trackResults(activeChnList).PRN];
    
    % These two lines help the skyPlot function. The satellites excluded
    % do to elevation mask will not "jump" to possition (0,0) in the sky
    % plot.
    navSolutions.channel.el(:, currMeasNr) = ...
        NaN(settings.numberOfChannels, 1);
    navSolutions.channel.az(:, currMeasNr) = ...
        NaN(settings.numberOfChannels, 1);
    
    %% Calculate the current sample number, corresponding satellites ======
    %% tranmitting time and raw receiver time =============================
    sampleNum=currMeasNr*settings.samplingFreq/settings.navSolRate...
        +settings.skipNumberOfSamples;
    transmitTime=...
        findTransTime(sampleNum,activeChnList,svTimeTable,trackResults);
    rxTime=max(transmitTime)+settings.startOffset/1000;
    
    %% Find pseudoranges ======================================================
    
    [navSolutions.channel.rawP(:, currMeasNr)] = calculatePseudoranges(...
        transmitTime,rxTime,activeChnList,settings);
    % old version of calculateP
%     navSolutions.channel.rawP(:, currMeasNr) = calculatePseudoranges(...
%             trackResults, ...
%             subFrameStart + 1000/settings.navSolRate * (currMeasNr-1), ...
%             activeChnList, settings);
    
    %% Find satellites positions and clocks corrections =======================
    
%     [satPositions, satClkCorr] = satpos(transmitTime, ...
%         [trackResults(activeChnList).PRN],eph);
        [satPositions, satClkCorr] = satpos(transmitTime(find(transmitTime>0)), ...
        [trackResults(activeChnList).PRN],eph);
    
    %% Find receiver position =================================================
    
    % 3D receiver position can be found only if signals from more than 3
    % satellites are available
    if length(activeChnList) > 3
%     if size(activeChnList, 2) > 3
        
        %=== Calculate receiver position ==================================
        freqforcal=zeros(1,length(activeChnList));
        for ii=1:length(activeChnList)
            freqforcal(ii)=trackResults(1,activeChnList(ii)).carrFreq(currMeasNr*1000/settings.navSolRate);
        end
        
        [xyzdt, ...
            navSolutions.channel.el(activeChnList, currMeasNr), ...
            navSolutions.channel.az(activeChnList, currMeasNr), ...
            navSolutions.DOP(:, currMeasNr)] = ...
            leastSquarePos(satPositions, ...
            navSolutions.channel.rawP(activeChnList, currMeasNr)' + satClkCorr * settings.c, ...
            freqforcal,settings);
        
%         [xyzdt, ...
%             navSolutions.channel.el(activeChnList, currMeasNr), ...
%             navSolutions.channel.az(activeChnList, currMeasNr), ...
%             navSolutions.DOP(:, currMeasNr)] = ...
%             leastSquarePos(satPositions, ...
%             navSolutions.channel.rawP(activeChnList, currMeasNr)' + satClkCorr * settings.c, ...
%             settings);
        
        %--- Save results -------------------------------------------------
        navSolutions.X(currMeasNr)  = xyzdt(1);
        navSolutions.Y(currMeasNr)  = xyzdt(2);
        navSolutions.Z(currMeasNr)  = xyzdt(3);
        navSolutions.dt(currMeasNr) = xyzdt(4);
        
        navSolutions.Vx(currMeasNr)  = xyzdt(5);
        navSolutions.Vy(currMeasNr)  = xyzdt(6);
        navSolutions.Vz(currMeasNr)  = xyzdt(7);
        navSolutions.ddt(currMeasNr) = xyzdt(8);
        
        % Update the satellites elevations vector
        satElev = navSolutions.channel.el(:, currMeasNr);
        
        %=== Correct pseudorange measurements for clocks errors ===========
        navSolutions.channel.correctedP(activeChnList, currMeasNr) = ...
            navSolutions.channel.rawP(activeChnList, currMeasNr) + ...
            satClkCorr' * settings.c - navSolutions.dt(currMeasNr);
        
        %% Coordinate conversion ==================================================
        
        %=== Convert to geodetic coordinates ==============================
        [navSolutions.latitude(currMeasNr), ...
            navSolutions.longitude(currMeasNr), ...
            navSolutions.height(currMeasNr)] = cart2geo(...
            navSolutions.X(currMeasNr), ...
            navSolutions.Y(currMeasNr), ...
            navSolutions.Z(currMeasNr), ...
            5);
        
        %=== Convert to UTM coordinate system =============================
        navSolutions.utmZone = findUtmZone(navSolutions.latitude(currMeasNr), ...
            navSolutions.longitude(currMeasNr));
        
        [navSolutions.E(currMeasNr), ...
            navSolutions.N(currMeasNr), ...
            navSolutions.U(currMeasNr)] = cart2utm(xyzdt(1), xyzdt(2), ...
            xyzdt(3), ...
            navSolutions.utmZone);
        
        % Compute the corrected receiver time
        navSolutions.rxTime(currMeasNr)=rxTime-navSolutions.dt(currMeasNr)/settings.c;
        % Record the sample number and raw receiver time
        navSolutions.absoluteSample(currMeasNr) =sampleNum;
        navSolutions.rawRxTime(currMeasNr)=rxTime;

        %DMA add - get the precise time of the first sample and the avg
        %clock rate for the file (skip first 5 samples and should be
        %enough to get over transients)
        if (currMeasNr == fix((min(lastSample) - ...
                settings.samplingFreq/settings.navSolRate-...
                settings.skipNumberOfSamples) /navStep))
            dmaTime=polyfit(navSolutions.absoluteSample(5:end)-settings.skipNumberOfSamples,navSolutions.rxTime(5:end),1);
            navSolutions.avgClock = 1/dmaTime(1);
            navSolutions.firstSampleTime = 1 * dmaTime(1) + dmaTime(2);
        end
    else % if size(activeChnList, 2) > 3
        %--- There are not enough satellites to find 3D position ----------
        disp(['   Measurement No. ', num2str(currMeasNr), ...
            ': Not enough information for position solution.']);

        %--- Set the missing solutions to NaN. These results will be
        %excluded automatically in all plots. For DOP it is easier to use
        %zeros. NaN values might need to be excluded from results in some
        %of further processing to obtain correct results.
        navSolutions.X(currMeasNr)           = NaN;
        navSolutions.Y(currMeasNr)           = NaN;
        navSolutions.Z(currMeasNr)           = NaN;
        navSolutions.dt(currMeasNr)          = NaN;
        navSolutions.DOP(:, currMeasNr)      = zeros(5, 1);
        navSolutions.latitude(currMeasNr)    = NaN;
        navSolutions.longitude(currMeasNr)   = NaN;
        navSolutions.height(currMeasNr)      = NaN;
        navSolutions.E(currMeasNr)           = NaN;
        navSolutions.N(currMeasNr)           = NaN;
        navSolutions.U(currMeasNr)           = NaN;
        navSolutions.rawRxTime               = NaN;
        navSolutions.absoluteSample          = NaN;
        navSolutions.rxTime                  = NaN;
        
        navSolutions.channel.az(activeChnList, currMeasNr) = ...
            NaN(1, length(activeChnList));
        navSolutions.channel.el(activeChnList, currMeasNr) = ...
            NaN(1, length(activeChnList));
        
        % TODO: Know issue. Satellite positions are not updated if the
        % satellites are excluded do to elevation mask. Therefore rasing
        % satellites will be not included even if they will be above
        % elevation mask at some point. This would be a good place to
        % update positions of the excluded satellites.
        
        
        
    end % if size(activeChnList, 2) > 3
    
end
toc
