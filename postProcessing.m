% Script postProcessing.m processes the raw signal from the specified data
% file (in settings) operating on blocks of 37 seconds of data.
%
% First it runs acquisition code identifying the satellites in the file,
% then the code and carrier for each of the satellites are tracked, storing
% the 1msec accumulations.  After processing all satellites in the 37 sec
% data block, then postNavigation is called. It calculates pseudoranges
% and attempts a position solutions. At the end plots are made for that
% block of data.


%--------------------------------------------------------------------------

%                         THE SCRIPT "RECIPE"
%
% The purpose of this script is to combine all parts of the software
% receiver.
%
% 1.1) Open the data file for the processing and seek to desired point.
%
% 2.1) Acquire satellites
%
% 3.1) Initialize channels (preRun.m).
% 3.2) Pass the channel structure and the file identifier to the tracking
% function. It will read and process the data. The tracking results are
% stored in the trackResults structure. The results can be accessed this
% way (the results are stored each millisecond):
% trackResults(channelNumber).XXX(fromMillisecond : toMillisecond), where
% XXX is a field name of the result (e.g. I_P, codePhase etc.)
%
% 4) Pass tracking results to the navigation solution function. It will
% decode navigation messages, find satellite positions, measure
% pseudoranges and find receiver position.
%
% 5) Plot the results.


%% Initialization =========================================================
disp ('Starting processing...');

[fid, message] = fopen(settings.fileName, 'rb');

%Initialize the multiplier to adjust for the data type
if (settings.fileType==1)
    dataAdaptCoeff=1;
else
    dataAdaptCoeff=2;
end

%If success, then process the data
if (fid > 0)

    % Move the starting point of processing. Can be used to start the
    % signal processing at any point in the data record (e.g. good for long
    % records or for signal processing in blocks).
    fseek(fid, dataAdaptCoeff*settings.skipNumberOfSamples, 'bof');

    %% Acquisition ============================================================

    % Do acquisition if it is not disabled in settings or if the variable
    % acqResults does not exist.
    if ((settings.skipAcquisition == 0))% || ~exist('acqResults', 'var'))

        % Find number of samples per spreading code
        samplesPerCode = round(settings.samplingFreq / ...
            (settings.codeFreqBasis / settings.codeLength));



        %--- Do the acquisition -------------------------------------------
        disp ('   Acquiring satellites...');

        % Read the required amount of data depending on the data file type
        % and the number of code period of coherent and non-coherent 
        % integration and invoke the acquisition function

        data = fread(fid, ...
            dataAdaptCoeff*(settings.acquisition.cohCodePeriods*settings.acquisition.nonCohSums+1)*samplesPerCode*6, ...
            settings.dataType)';

        if (dataAdaptCoeff==2)
            data1=data(1:2:end);
            data2=data(2:2:end);
            data=data1 + 1i .* data2;
        end

        acqResults = acquisition(data, settings);


        % Plot the acquisition results
        plotAcquisition(acqResults);
    elseif(settings.skipAcquisition == 1)
            disp('skipAcquistion==1');
            load acqresults
            plotAcquisition(acqResults);
    end

    %% Initialize channels and prepare for the run ============================

    % Start further processing only if a GNSS signal was acquired (the
    % field FREQUENCY will be set to 0 for all not acquired signals)
    if (any(acqResults.peakMetric>settings.acqThreshold))
        channel = preRun(acqResults, settings);
        showChannelStatus(channel, settings);
    else
        % No satellites to track, exit
        disp('No GNSS signals detected, signal processing finished.');
        trackResults = [];
        return;
    end

    %% Track the signal =======================================================
    if ~settings.VLLen
       
        startTime = now;
        disp (['   Tracking started at ', datestr(startTime)]);

        % Process all channels for given data block
       
        [trackResults, channel] = tracking(fid, channel, settings);
       
        % Close the data file
        fclose(fid);

        disp(['   Tracking is over (elapsed time ', ...
            datestr(now - startTime, 13), ')'])

        % Auto save the acquisition & tracking results to a file to allow
        % running the positioning solution afterwards.
        disp('   Saving Acq & Tracking results to file "trackingResults.mat"')
        save('trackingResults', ...
            'trackResults', 'settings', 'acqResults', 'channel');
%%%%%%%%%% when delete, uncomment above
    load trackingResults

    else
        disp('skip scalar tracking, load exisiting results')

        load 'trackingResults'
        settings.VLLen = 1;
    end

    %% Calculate navigation solutions =========================================
    
    if ~settings.VLLen
        disp('   Calculating navigation solutions...');
        [navSolutions, eph, svTimeTable,activeChnList] = postNavigation(trackResults, settings);
        save('navSolutions','navSolutions','eph','svTimeTable','activeChnList')
        load 'navSolutions'
    else
        disp('skip postNavigation, load existing navSolution')
        %input the path and file name of the existing navSolution saved
        %from scalar tracking and calculation
        load 'navSolutions'
        disp('start vector tracking')
        return;
%         load 'C:\Users\gps\Desktop\SihaoZ\dynamic\9_13_2010_14h_53min_25s\navSolutionscarfig8allfading190spll10dll1'
                %start vector tracking
        [trackResults_v, channel] = trackingv(fid, channel,trackResults,navSolutions,eph,activeChnList,svTimeTable, settings);
        save('trackResults_v','trackResults_v');
    end

    disp('   Processing is complete for this data block');
    %% Plot all results ===================================================
    disp ('   Ploting results...');
    if settings.plotTracking
        plotTracking(1:settings.numberOfChannels, trackResults, settings);
    end

    plotNavigation(navSolutions, settings);

    disp('Post processing of the signal is over.');

else
    % Error while opening the data file.
    error('Unable to read file %s: %s.', settings.fileName, message);
end % if (fid > 0)