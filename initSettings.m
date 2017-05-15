function settings = initSettings()
%Functions initializes and saves settings. Settings can be edited inside of
%the function, updated from the command line or updated using a dedicated
%GUI - "setSettings".
%
%All settings are described inside function code.
%
%settings = initSettings()
%
%   Inputs: none
%
%   Outputs:
%       settings     - Receiver settings (a structure).

%--------------------------------------------------------------------------

%% Processing settings ====================================================
% Number of milliseconds to be processed used 36000 + any transients (see
% below - in Nav parameters) to ensure nav subframes are provided
settings.msToProcess        = 36000;        %[ms]

% Number of channels to be used for signal processing
settings.numberOfChannels   = 8;

% Move the starting point of processing. Can be used to start the signal
% processing at any point in the data record (e.g. for long records). fseek
% function is used to move the file read point, therefore advance is byte
% based only. For Real sample files it skips the number of bytes as indicated
% here. For I/Q files it skips twice the number of bytes as indicated here
% to consider both I and Q samples
settings.skipNumberOfSamples     = 0*2.5e6;
settings.skipNumberOfSamples     = 1*53e6;
settings.skipNumberOfBytes     = 0;

%% Raw signal file name and other parameter ===============================
% This is a "default" name of the data file (signal record) to be used in
% the post-processing mode

% settings.fileName           = ...
%     'D:\forZhu\8_30_2010_15h_15min_48s\8_30_2010_15h_15min_48s_SDR.bin'
%     %'C:\Users\gps\Desktop\SihaoZ\dynamic\9_13_2010_14h_53min_25s\9_13_2010_14h_53min_25s_SDR.bin';
%     %'C:\Nordnav-Rx_v3-6-1\h22m49up.sim';
% settings.fileName           =  'C:\gnss0.bin';
% settings.fileName = 'C:\Users\AltBOC\Documents\MATLAB\sdr\data\201404150757amUTCG\v2\gnsa14.bin';
% settings.fileName = 'C:\Users\AltBOC\Documents\MATLAB\gnsa14.bin';
% settings.fileName = 'D:\000Unicorn\IFRecords\GN3Sv3\dump_mode1.bin'; %
% 60s, 20170508, 0:46, UTC+1, BST %% I & Q samples
% settings.fileName = 'D:\000Unicorn\IFRecords\GN3Sv3\dump_mode3.bin'; %
% 60s, 20170508, 2:20, UTC+1, BST %% only real-valued samples
settings.fileName = 'D:\000Unicorn\IFRecords\NUT4NT\dump4ch_ch2.bin';
settings.fileName = 'I:\000_NUT4NT\GPS_L1_2016_10_22_140330.bin';
    
    
% Data type used to store one sample
settings.dataType           = 'schar';

% File Types
%1 - 8 bit real samples S0,S1,S2,...
%2 - 8 bit I/Q samples I0,Q0,I1,Q1,I2,Q2,...
settings.fileType           = 1;

% Intermediate, sampling and code frequencies

% settings.IF                 = 0e6;%-1e6;%10;      %[Hz]
% settings.samplingFreq       = 2.5e6;%5e6;%2.048e6;    %[Hz]
% % % % settings.IF                 = 38400;        %[Hz]
% % % % settings.samplingFreq       = 16.3676e6/2;  %[Hz]
% % % % settings.codeFreqBasis      = 1.023e6;      %[Hz]

%% mode 3
settings.IF                 = 14.58e6;        %[Hz]
settings.samplingFreq       = 53e6;  %[Hz]
settings.codeFreqBasis      = 1.023e6;      %[Hz]


% Define number of chips in a code period
settings.codeLength         = 1023;

%% Acquisition settings ===================================================
% Skips acquisition in the script postProcessing.m if set to 1
settings.skipAcquisition    = 0;
% List of satellites to look for. Some satellites can be excluded to speed
% up acquisition
% % % % settings.acqSatelliteList   = [1:29 31 32];%[1:32];         %[PRN numbers]PRN32 was just launched (as of 2014 Apr 15) so it was not set active yet:accuracy indicated by decoded message:11
settings.acqSatelliteList   = 1:32;
% Band around IF to search for satellite signal. Depends on max Doppler
settings.acqSearchBand      = 14;           %[kHz]
% Threshold for the signal presence decision rule
settings.acqThreshold       = 2.5;
% No. of code periods for coherent integration (multiple of 2)
settings.acquisition.cohCodePeriods=2;%10;%
% No. of non-coherent summations
settings.acquisition.nonCohSums=4;

%% Tracking loops settings ================================================
settings.enableFastTracking     = 0;

% Code tracking loop parameters
settings.dllDampingRatio         = 0.7;
settings.dllNoiseBandwidth       = 1;       %[Hz]
settings.dllCorrelatorSpacing    = 0.5;     %[chips]
% Carrier tracking loop parameters
settings.pllDampingRatio         = 0.7;
settings.pllNoiseBandwidth       = 6.5;      %[Hz]
settings.fllDampingRatio         = 0.7;
settings.fllNoiseBandwidth       = 10;      %[Hz]



%% Navigation solution settings ===========================================

% Rate for calculating pseudorange and position
settings.navSolRate         = 10;            %[Hz]
settings.navSolPeriod       = 1000/settings.navSolRate; %ms

% Elevation mask to exclude signals from satellites at low elevation
settings.elevationMask      = 10;           %[degrees 0 - 90]
% Enable/dissable use of tropospheric correction
settings.useTropCorr        = 0;            % 0 - Off
% 1 - On

% True position of the antenna in UTM system (if known). Otherwise enter
% all NaN's and mean position will be used as a reference .
settings.truePosition.E     = nan;
settings.truePosition.N     = nan;
settings.truePosition.U     = nan;

%% Plot settings ==========================================================
% Enable/disable plotting of the tracking results for each channel
settings.plotTracking       = 1;
% 0 - Off
% 1 - On

%% Constants ==============================================================
settings.c                  = 299792458;    % The speed of light, [m/s]
settings.startOffset        = 68.802;       %[ms] Initial sign. travel time

%% CNo Settings============================================================
% Accumulation interval in Tracking (in Sec)
settings.CNo.accTime=0.001;
% Show C/No during Tracking;1-on;0-off;
% settings.CNo.enableVSM=1;
settings.CNo.enableVSM=0;
% Accumulation interval for computing VSM C/No (in ms)
settings.CNo.VSMinterval=400;
% Accumulation interval for computing PRM C/No (in ms)
settings.CNo.PRM_K=200;
% No. of samples to calculate narrowband power;
% Possible Values for M=[1,2,4,5,10,20];
% K should be an integral multiple of M i.e. K=nM
settings.CNo.PRM_M=20;
% Accumulation interval for computing MOM C/No (in ms)
settings.CNo.MOMinterval=200;
% Enable/disable the C/No plots for all the channels
% 0 - Off ; 1 - On;
settings.CNo.Plot = 1;
%Enable vector tracking when 1, otherwise scalar tracking.
% settings.VLLen = 1;
settings.VLLen = 0;