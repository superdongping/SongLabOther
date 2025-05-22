clc; clear; close all;

% mouse_colony_analysis.m
% Analyze mouse colony: summary stats, combined subplots, and pie charts.

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

% Loop through each line
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
            if strcmpi(T.(dispCol){j},'Mating') && monthsAge>6
                T.Note(j) = "might need to update the pairing cage";
            end
        end
    end
    
    % Subplot and histogram
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
writetable(T,'MouseColony_Analyzed.xlsx');

% Display summary
disp('Summary of cages and mice per mouse line:');
disp(summaryData);
fprintf('Analysis complete. Results saved to MouseColony_Analyzed.xlsx\n');

% --- Additional Figures ---
% Pie chart: number of mice per mouse line
figure('Name','Mice Count per Line');
labels = num2str(summaryData.NumMice);
pie(summaryData.NumMice, labels);
legend(summaryData.MouseLine,'Location','eastoutside');
title('Mice Count by Mouse Line');

% Pie chart: number of cages per mouse line
figure('Name','Cage Count per Line');
labels = num2str(summaryData.NumCages);
pie(summaryData.NumCages,labels);
legend(summaryData.MouseLine,'Location','eastoutside');
title('Cage Count by Mouse Line');
