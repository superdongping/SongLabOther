clc;
close all;
clear all;
%% Load Fiber Photometry CSV Data
[file, path] = uigetfile('*.csv', 'Select the fiber photometry CSV file');
if isequal(file, 0)
    disp('User cancelled file selection.');
    return;
end
csvFile = fullfile(path, file);

% Read the CSV file into a table. Import time as strings.
data = readtable(csvFile, 'TextType', 'string');

% Extract the time column and convert it to a duration object.
% The expected format is d:hh:mm:ss.FFFFFFF.
TimeStamp = data{:, 1};

% Convert duration to seconds using the provided snippet
if isduration(TimeStamp)
    timeVector = seconds(TimeStamp) - seconds(TimeStamp(1));
else
    timeVector = seconds(TimeStamp - TimeStamp(1));
end
timeVector = timeVector(:);  % Ensure it's a column vector

% Convert time to hours for plotting on the x-axis
time_hours = timeVector / 3600;

% Extract the LED signals (assumed: column 2 is 410 and column 3 is 470)
LED410 = data{:,2};
LED470 = data{:,3};

% Compute the ratio (470/410)
ratio = LED470 ./ LED410;

%% Plot Raw Data (Normalized)
% Normalize raw data by subtracting its minimum value
raw_normalized = ratio - min(ratio);

% Downsample data if needed
maxPoints = 10000;
numPoints = numel(time_hours);
if numPoints > maxPoints
    ds_factor = ceil(numPoints / maxPoints);
    time_hours_ds = time_hours(1:ds_factor:end);
    raw_normalized_ds = raw_normalized(1:ds_factor:end);
else
    time_hours_ds = time_hours;
    raw_normalized_ds = raw_normalized;
end

% Plot raw normalized data
figure;
plot(time_hours_ds, raw_normalized_ds, '-b');
xlabel('Time (hours)');
ylabel('Normalized Ratio (470/410)');
title('Raw Normalized Fiber Photometry Signal over 24 Hours');
grid on;
xlim([0 24]);
xticks(0:1:24);

%% Apply Baseline Correction and Plot Corrected Data
% Apply baseline correction using a linear detrending (see function below)
ratio_corrected = applyBaselineCorrection(ratio, timeVector);

% Normalize the baseline corrected data by subtracting its minimum value
corrected_normalized = ratio_corrected - min(ratio_corrected);

% Downsample baseline corrected data (using same factor as before)
if numPoints > maxPoints
    corrected_normalized_ds = corrected_normalized(1:ds_factor:end);
else
    corrected_normalized_ds = corrected_normalized;
end

%% Compute the Slow Component (Baseline)
% Here, we extract the slow component using a cubic fit to the original ratio signal.
p = polyfit(timeVector, ratio, 3);
baseline = polyval(p, timeVector);
% Normalize the baseline (slow component) by subtracting its minimum value
baseline_norm = baseline - min(baseline);

% Downsample the slow component for plotting
if numPoints > maxPoints
    baseline_norm_ds = baseline_norm(1:ds_factor:end);
else
    baseline_norm_ds = baseline_norm;
end

%% Plot Baseline Corrected Data with Overlaid Slow Component
figure;
plot(time_hours_ds, corrected_normalized_ds, '-r', 'DisplayName', 'Baseline Corrected');
hold on;
plot(time_hours_ds, baseline_norm_ds, '-k', 'LineWidth', 2, 'DisplayName', 'Slow Component');
xlabel('Time (hours)');
ylabel('Normalized Signal (470/410)');
title('Baseline Corrected Signal with Slow Component over 24 Hours');
legend('show');
grid on;
xlim([0 24]);
xticks(0:1:24);

%% Compute Hourly Summaries (AUC and Average Value) from Baseline Corrected Data
% Define hourly bins (0 to 24 hours)
hour_edges = 0:1:24;
nBins = length(hour_edges) - 1;
AUC = zeros(nBins, 1);
avg_val = zeros(nBins, 1);

for i = 1:nBins
    % Find indices for the current hour bin
    bin_idx = time_hours >= hour_edges(i) & time_hours < hour_edges(i+1);
    if any(bin_idx)
        % Compute average normalized value for this hour
        avg_val(i) = mean(corrected_normalized(bin_idx));
        % Compute area under the curve using numerical integration (trapz)
        AUC(i) = trapz(time_hours(bin_idx), corrected_normalized(bin_idx));
    else
        avg_val(i) = NaN;
        AUC(i) = NaN;
    end
end

%% Plot Hourly Summary in a New Figure
figure;
subplot(2,1,1);
bar(hour_edges(1:end-1)+0.5, AUC, 'FaceColor', [0.2, 0.6, 0.8]);
xlabel('Hour');
ylabel('AUC');
title('Area Under the Curve per Hour');
xlim([0 24]);
xticks(0:1:24);

subplot(2,1,2);
bar(hour_edges(1:end-1)+0.5, avg_val, 'FaceColor', [0.8, 0.4, 0.4]);
xlabel('Hour');
ylabel('Average Normalized Signal');
title('Average Value per Hour');
xlim([0 24]);
xticks(0:1:24);

%% Save Hourly Summary to Excel File
% Create a table with Hour, AUC, and Average Value
hour_centers = (hour_edges(1:end-1) + 0.5)';
summaryTable = table(hour_centers, AUC, avg_val, ...
    'VariableNames', {'Hour', 'AUC', 'Average_Normalized_Signal'});

% Write the table to an Excel file named "NPY_24H_summary.xlsx" in the same folder as the CSV file
excelFile = fullfile(path, 'NPY_24H_summary.xlsx');
writetable(summaryTable, excelFile);
fprintf('Hourly summary saved to: %s\n', excelFile);

%% Baseline Correction Function
function correctedSignal = applyBaselineCorrection(signal, timeVector)
% This function removes a slow drift (baseline) from the signal using linear detrending.
% Fit a linear polynomial trend to the signal.
validIdx = isfinite(timeVector) & isfinite(signal);
if any(~validIdx)
    warning('Some data points are non-finite and will be omitted.');
end
p = polyfit(timeVector(validIdx), signal(validIdx), 1);


% p = polyfit(timeVector, signal, 1);
% Calculate the baseline from the fit.
baseline = polyval(p, timeVector);
% Subtract the baseline from the original signal.
correctedSignal = signal - baseline;
end
