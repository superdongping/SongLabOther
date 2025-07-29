%% Batch Processing of ABF Files for Gap-Free EPSC Analysis with AUC Filtering
% This script allows the user to select a folder containing ABF files,
% processes each file to detect EPSC events across the entire recording,
% and applies an AUC-based criterion to exclude noise.

clc; clear; close all;

%% Let the user choose a folder containing ABF files
dataFolder = uigetdir(pwd, 'Select folder containing ABF files');
if isequal(dataFolder, 0)
    disp('No folder selected. Exiting.');
    return;
end

%% List all .abf files in the selected folder
abfFiles = dir(fullfile(dataFolder, '*.abf'));
if isempty(abfFiles)
    error('No ABF files found in the selected folder.');
end

%% Ask user if they want to display detection plots
choice = questdlg('Show detection plots for each sweep?','EPSC Detection','Yes','No','No');
showPlots = strcmp(choice,'Yes');

%% Analysis parameters (adjust as needed)
fs = 10000;                    % Sampling rate (Hz)
dt = 1/fs;                     % Sampling interval (s)

% Filter parameters
cutoff = 100;                  % Low-pass filter cutoff (Hz)
order = 4;                     % Butterworth filter order

% Baseline and smoothing
polyDegree = 4;                % Degree for polynomial baseline fit
smoothingWindow = 50;          % Smoothing window (samples)

% Peak detection
minEventAmplitude = 5;        % Min peak amplitude (pA)
minPeakDistanceSamples = round(0.1 * fs);  % Min peak-to-peak interval (samples)

% AUC filtering parameters
preWindow_s = 0.005;           % 5 ms before peak
postWindow_s = 0.01;           % 10 ms after peak
preWindow = round(preWindow_s * fs);
postWindow = round(postWindow_s * fs);
minEventAUC = 0.01;               % Minimum area under curve (pA*s) to include event

% Design Butterworth filter
[b, a] = butter(order, cutoff/(fs/2), 'low');

%% Initialize master summary across all files
masterSummary = table();

%% Loop over each ABF file
for fIdx = 1:numel(abfFiles)
    abfName = fullfile(abfFiles(fIdx).folder, abfFiles(fIdx).name);
    [~, baseName, ~] = fileparts(abfName);
    disp(['Processing file: ', abfName]);

    % Load ABF data (requires abfload on path)
    [data, ~, ~] = abfload(abfName);
    [nSamples, ~, nSweeps] = size(data);
    timeVec = (0:nSamples-1) * dt;  % time vector (row)

    % Per-file storage
    summaryData = zeros(nSweeps, 5);  % [Sweep, Count, MeanAmp, MeanISI, MeanAUC]
    fileRawEPSC = table();
    firstSweepCorr = [];
    firstSweepValidPeaks = [];

    %% Loop through sweeps
    for sw = 1:nSweeps
        trace = data(:,1,sw);

        % Filtering and baseline correction
        filtTrace = filtfilt(b, a, trace);
        tCol = timeVec.';
        p = polyfit(tCol, filtTrace, polyDegree);
        baseline = polyval(p, tCol);
        corrTrace = filtTrace - baseline;
        smoothTrace = smoothdata(corrTrace, 'movmean', smoothingWindow);

        % Detect candidate EPSC peaks
        [peakVals, peakLocs] = findpeaks(-smoothTrace, ...
            'MinPeakHeight', minEventAmplitude, ...
            'MinPeakDistance', minPeakDistanceSamples);
        amplitudes = -peakVals;

        % Compute AUC around each candidate
        nPeaks = numel(peakLocs);
        aucValsAll = nan(nPeaks,1);
        for pi = 1:nPeaks
            idx0 = peakLocs(pi);
            idxStart = max(1, idx0 - preWindow);
            idxEnd = min(nSamples, idx0 + postWindow);
            tWin = timeVec(idxStart:idxEnd);
            xWin = corrTrace(idxStart:idxEnd);
            aucValsAll(pi) = trapz(tWin, abs(xWin));
        end

        % Filter peaks by AUC minimum
        validIdx = aucValsAll >= minEventAUC;
        peakLocs = peakLocs(validIdx);
        amplitudes = amplitudes(validIdx);
        aucVals = aucValsAll(validIdx);

        % Compute event times and ISIs
        eventTimes = timeVec(peakLocs).';
        if numel(eventTimes) > 1
            isi = [NaN; diff(eventTimes)];
        else
            isi = NaN;
        end

        % Store first-sweep valid peaks for PNG
        if sw == 1
            firstSweepCorr = corrTrace;
            firstSweepValidPeaks = peakLocs;
        end

        % Summary metrics including AUC
        count = numel(amplitudes);
        if count > 0
            meanAmp = mean(amplitudes, 'omitnan');
        else
            meanAmp = NaN;
        end
        if numel(eventTimes) > 1
            meanISI = mean(diff(eventTimes));
        else
            meanISI = NaN;
        end
        if ~isempty(aucVals)
            meanAUC = mean(aucVals, 'omitnan');
        else
            meanAUC = NaN;
        end
        summaryData(sw,:) = [sw, count, meanAmp, meanISI, meanAUC];

        % Append raw EPSC events with AUC
        if count > 0
            nEvt = numel(amplitudes);
            temp = table(repmat(sw,nEvt,1), (1:nEvt).', eventTimes, amplitudes, isi, aucVals, ...
                'VariableNames', {'Sweep','EventNumber','Time_s','Amplitude_pA','ISI_s','AUC_pA_s'});
            fileRawEPSC = [fileRawEPSC; temp];
        end

        % On-screen plot if requested
        if showPlots
            figure('Name', sprintf('%s - Sweep %d', baseName, sw));
            plot(timeVec, corrTrace, 'b'); hold on;
            plot(timeVec(firstSweepValidPeaks), corrTrace(firstSweepValidPeaks), 'r*');
            xlabel('Time (s)'); ylabel('Corrected Current (pA)');
            title(sprintf('%s - Sweep %d EPSC Detection (AUC filtered)', baseName, sw));
            grid on; pause(0.5);
        end
    end

    %% Create and save per-file summary and raw tables
    summaryTable = array2table(summaryData, ...
        'VariableNames', {'Sweep','EPSC_Count','EPSC_MeanAmp_pA','EPSC_MeanISI_s','EPSC_MeanAUC_pA_s'});
    outName = fullfile(dataFolder, [baseName, '_EPSC_summary.xlsx']);
    writetable(summaryTable, outName, 'Sheet', 1);
    writetable(fileRawEPSC, outName, 'Sheet', 'Raw_EPSC');
    disp(['Saved summary and raw events: ', outName]);

    %% Save detection overlay PNG for first sweep
    if ~isempty(firstSweepCorr)
        fig = figure('Visible','off');
        plot(timeVec, firstSweepCorr, 'b'); hold on;
        peakTimes = timeVec(firstSweepValidPeaks);
        plot(peakTimes, firstSweepCorr(firstSweepValidPeaks), 'r*');
        xlabel('Time (s)'); ylabel('Corrected Current (pA)');
        title(sprintf('%s - Sweep 1 EPSC Detection (AUC filtered)', baseName));
        grid on;
        pngName = fullfile(dataFolder, [baseName, '_EPSC_detection.png']);
        saveas(fig, pngName);
        close(fig);
        disp(['Saved detection PNG: ', pngName]);
    end

    %% Append to master summary
    totalEvents = sum(summaryTable.EPSC_Count);
    avgAmp = mean(summaryTable.EPSC_MeanAmp_pA, 'omitnan');
    avgISI = mean(summaryTable.EPSC_MeanISI_s, 'omitnan');
    avgAUC = mean(summaryTable.EPSC_MeanAUC_pA_s, 'omitnan');
    masterSummary = [masterSummary; table({abfName}, totalEvents, avgAmp, avgISI, avgAUC, ...
        'VariableNames', {'FileName','TotalEPSC','AvgAmp_pA','AvgISI_s','AvgAUC_pA_s'})];
end

%% Save master summary across all files
masterOut = fullfile(dataFolder, 'Master_EPSC_Summary.xlsx');
writetable(masterSummary, masterOut, 'Sheet', 1);
disp(['Master summary saved: ', masterOut]);
