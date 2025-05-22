clc; clear; close all;

% mouse_colony_analysis.m
% Analyze mouse colony: summary stats, combined subplots, and detailed pie charts.

% Define file name
filename = 'SoftMouse.NET-Cage List-PingDong2025-05-22 0939.xlsx';

% Read table preserving original names
opts = detectImportOptions(filename, 'VariableNamingRule','preserve');
T = readtable(filename, opts);
varNames = T.Properties.VariableNames;

% Identify key columns
lineCol    = varNames{contains(varNames,'Mouseline','IgnoreCase',true)};   % Mouse line
tagCol     = varNames{contains(varNames,'Tags','IgnoreCase',true)};        % Tags with DOB
dispCol    = varNames{contains(varNames,'Disposition','IgnoreCase',true)};  % Cage disposition
numMCols   = varNames(contains(varNames,'#') & contains(varNames,'Mice','IgnoreCase',true));
numMiceCol = numMCols{1};                                                % Number of mice

% Convert number-of-mice to numeric if needed
rawNumMice = T.(numMiceCol);
if iscell(rawNumMice) || isstring(rawNumMice)
    numMiceData = str2double(rawNumMice);
else
    numMiceData = rawNumMice;
end

% Date thresholds
todayDate    = datetime('today');
sixMonthsAgo = todayDate - calmonths(6);

% Initialize Note column
T.Note = repmat(string(missing), height(T),1);

% Unique lines and summary
enumLines   = unique(string(T.(lineCol)));
nLines      = numel(enumLines);
summaryData = table('Size',[nLines,3], 'VariableTypes',{'string','double','double'}, 'VariableNames',{'MouseLine','NumCages','NumMice'});

% Prepare subplot layout: square-ish
gridRows = ceil(sqrt(nLines));
gridCols = ceil(nLines/gridRows);
figure('Name','Age Distributions','Units','normalized','Position',[0.1 0.1 0.8 0.8]);

% Loop through each line: histogram and summary
for i = 1:nLines
    line    = enumLines(i);
    idxLine = string(T.(lineCol))==line;
    % Summary counts
    summaryData.MouseLine(i) = line;
    summaryData.NumCages(i)  = sum(idxLine);
    summaryData.NumMice(i)   = sum(numMiceData(idxLine));
    
    % Collect ages in months
    ageMonth = [];
    for j = find(idxLine)'
        tokens = regexp(T.(tagCol){j},'\d{2}-\d{2}-\d{4}','match');
        for k = 1:numel(tokens)
            dob       = datetime(tokens{k},'InputFormat','MM-dd-yyyy');
            monthsAge = days(todayDate-dob)/30;
            ageMonth(end+1) = monthsAge; %#ok<SAGROW>
            % Flag old mating cage
            if strcmpi(T.(dispCol){j},'Mating') && monthsAge > 6
                T.Note(j) = "might need to update the pairing cage";
            end
        end
    end
    
    % Subplot: age distribution histogram
    subplot(gridRows,gridCols,i);
    if ~isempty(ageMonth)
        maxM = ceil(max(ageMonth));
        edges = 0:1:maxM;
        histogram(ageMonth,edges);
        xticks(1:maxM);
        xticklabels(arrayfun(@(m) sprintf('%d month',m),1:maxM,'UniformOutput',false));
        xtickangle(45);
    end
    title(char(line),'Interpreter','none');
    xlabel('Age (months)');
    ylabel('Count');
    grid on;
end

% Write updated table
tblOut = 'MouseColony_Analyzed.xlsx';
writetable(T,tblOut);

% Display summary
disp('Summary of cages and mice per mouse line:');
disp(summaryData);
fprintf('Analysis complete. Results saved to %s\n',tblOut);

% --- Combined Pie Chart: Mice & Cages by Line & Disposition ---
% Compute mating and stock counts per line
matingMice = zeros(1,nLines);
stockMice  = zeros(1,nLines);
matingCages = zeros(1,nLines);
stockCages  = zeros(1,nLines);
for i = 1:nLines
    idxLine = string(T.(lineCol))==enumLines(i);
    matingIdx = idxLine & strcmpi(T.(dispCol),'Mating');
    stockIdx  = idxLine & ~strcmpi(T.(dispCol),'Mating');
    matingMice(i)  = sum(numMiceData(matingIdx));
    stockMice(i)   = sum(numMiceData(stockIdx));
    matingCages(i) = sum(matingIdx);
    stockCages(i)  = sum(stockIdx);
end

% Prepare data and labels for pie
data = [matingMice, stockMice];
labels = cell(1,2*nLines);
for i = 1:nLines
    labels{2*i-1} = sprintf('%s Mating', char(enumLines(i)));
    labels{2*i}   = sprintf('%s Stock', char(enumLines(i)));
end

% Generate distinct base colors for lines
baseColors = jet(nLines);
% Build colormap: mating=base; stock=lightened
cmap = zeros(2*nLines,3);
for i = 1:nLines
    cmap(2*i-1,:) = baseColors(i,:);
    % light variant
    cmap(2*i,:) = baseColors(i,:) + (1-baseColors(i,:))*0.8;
end

% Plot combined pie for mice counts
titleText = 'Mice Distribution by Line and Disposition';
figure('Name',titleText);
pie(data);
colormap(cmap);
legend(labels,'Location','eastoutside');
title(titleText);

% Plot combined pie for cage counts
data2 = [matingCages, stockCages];
titleText2 = 'Cage Distribution by Line and Disposition';
figure('Name',titleText2);
pie(data2);
colormap(cmap);
legend(labels,'Location','eastoutside');
title(titleText2);
