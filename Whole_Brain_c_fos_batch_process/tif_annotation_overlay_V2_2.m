%% Batch Processing Script (Index-Based Matching of TIFs and ROI Folders)

clear; clc; close all;

%% 0. Folders
% Folder where TIFF files reside
imageFolder = pwd;
% Base folder containing subfolders with ROI zip files
roiBaseFolder = pwd;
% Path to Allen Brain Atlas JSON
jsonFilePath = 'Allen_Brain_Atlas_database.json';

if ~isfile(jsonFilePath)
    error('JSON file "%s" not found.', jsonFilePath);
end

%% 1. Gather the TIFF files
tifFiles = dir(fullfile(imageFolder, '*.tif'));
% Sort them by name (or by date if you prefer). This ensures consistent order.
[~, sortIdx] = sort({tifFiles.name});
tifFiles = tifFiles(sortIdx);

numTifs = numel(tifFiles);
fprintf('Found %d TIFF files in: %s\n', numTifs, imageFolder);

%% 2. Gather the ROI subfolders
% Get list of subfolders
roiFolders = dir(roiBaseFolder);
roiFolders = roiFolders([roiFolders.isdir]);  % keep only directories
roiFolders = roiFolders(~ismember({roiFolders.name}, {'.','..'})); % remove '.' and '..'

% Convert the folder names to numbers
folderNums = zeros(numel(roiFolders), 1);
for i = 1:numel(roiFolders)
    folderNums(i) = str2double(roiFolders(i).name); 
    % str2double() will give NaN if the name is not purely numeric
end

% Sort numerically by these folder numbers
[~, sortIdx] = sort(folderNums);
roiFolders = roiFolders(sortIdx);

numFolders = numel(roiFolders);
fprintf('Found %d ROI subfolders in: %s\n', numFolders, roiBaseFolder);

% Sanity check: we expect the same number of TIF files and subfolders
if numTifs ~= numFolders
    error('Mismatch: %d TIFF files but %d subfolders. They must match 1-to-1.', numTifs, numFolders);
end

%% 3. Load and parse the Allen Brain Atlas JSON
jsonStr = fileread(jsonFilePath);
jsonData = jsondecode(jsonStr);
% Create a map to store region info as structs with fields 'acronym' and 'name'
regionInfoMap = containers.Map('KeyType','double','ValueType','any');
for i = 1:length(jsonData.msg)
    regionInfoMap = traverseJson(jsonData.msg(i), regionInfoMap);
end

%% 4. Loop through each TIF file and its matching ROI folder by index
for iFile = 1:numTifs
    try
        % --- TIF file to process ---
        tifImageFile = fullfile(imageFolder, tifFiles(iFile).name);
        fprintf('\nProcessing TIF #%d: %s\n', iFile, tifFiles(iFile).name);

        % --- ROI folder that corresponds to this TIF ---
        thisRoiFolder = fullfile(roiBaseFolder, roiFolders(iFile).name);
        fprintf('Matching ROI folder: %s\n', thisRoiFolder);

        % Within that folder, assume there is a single zip file:
        roiZipPath = fullfile(thisRoiFolder, ...
            'ABBA-RoiSet-Adult Mouse Brain - Allen Brain Atlas V3p1.zip');
        if ~isfile(roiZipPath)
            error('Could not find ROI zip: %s', roiZipPath);
        end

        %% 4a. Read the 3-channel TIF
        info = imfinfo(tifImageFile);
        numPages = numel(info);
        if numPages < 3
            error('TIF image does not have >=3 channels: %s', tifFiles(iFile).name);
        end

        img = zeros(info(1).Height, info(1).Width, numPages, 'double');
        for k = 1:numPages
            img(:,:,k) = double(imread(tifImageFile, k, 'Info', info));
        end

        % Clip intensities to [0, 5000], then normalize
        DAPI = min(max(img(:,:,1), 0), 5000) / 5000;  % channel 1
        cfos = min(max(img(:,:,3), 0), 5000) / 5000;  % channel 3

        % Create merged RGB: c-fos in R/G, DAPI in B
        alpha = cfos;
        merged = cat(3, alpha, alpha, (1 - alpha).*DAPI + alpha);

        %% 4b. Unzip ROI files to a temporary location
        tempRoiDir = fullfile(tempdir, ['ROI_tmp_' roiFolders(iFile).name]);
        if ~exist(tempRoiDir, 'dir')
            mkdir(tempRoiDir);
        end
        % Clean out old contents if needed:
        delete(fullfile(tempRoiDir, '*.roi'));
        unzip(roiZipPath, tempRoiDir);

        % Gather .roi files
        roiFiles = dir(fullfile(tempRoiDir, '*.roi'));
        if isempty(roiFiles)
            error('No .roi files found in unzipped folder: %s', tempRoiDir);
        end

        %% 4c. Parse the ROI files
        roiData = struct('fileName', {}, 'numericCoords', {}, 'rawData', {}, 'method', {});
        for k = 1:length(roiFiles)
            roiFile = fullfile(tempRoiDir, roiFiles(k).name);
            tempData = readImageJROI(roiFile);  % <-- your custom function
            numericCoords = [];
            method = '';

            if isfield(tempData, 'mnCoordinates')
                coords = tempData.mnCoordinates;
                if ischar(coords)
                    numericCoords = sscanf(coords, '%f', [2, Inf])';
                else
                    numericCoords = coords;
                end
                if ~isempty(numericCoords) && ~isequal(numericCoords(1,:), numericCoords(end,:))
                    numericCoords(end+1,:) = numericCoords(1,:);
                end
                method = 'mnCoordinates';

            elseif isfield(tempData, 'vfShapes')
                vfShapesData = tempData.vfShapes;
                polygons = parseVfShapes(vfShapesData);
                numericCoords = polygons;
                method = 'vfShapes';

            elseif isfield(tempData, 'roiPosition')
                pos = tempData.roiPosition; % [x, y, width, height]
                numericCoords = [ pos(1), pos(2);
                                  pos(1)+pos(3), pos(2);
                                  pos(1)+pos(3), pos(2)+pos(4);
                                  pos(1), pos(2)+pos(4);
                                  pos(1), pos(2)];
                method = 'roiPosition';
            else
                fprintf('No recognized coords in ROI file: %s\n', roiFiles(k).name);
            end

            roiData(k).fileName = roiFiles(k).name;
            roiData(k).numericCoords = numericCoords;
            roiData(k).rawData = tempData;
            roiData(k).method = method;
        end

        %% 4d. c-fos Positive Detection and ROI Quantification
        raw_image = cfos;
        % Example segmentation pipeline
        threshold = 0.5;
        I_BW = imbinarize(raw_image, threshold);
        I_BW_m = medfilt2(I_BW, [3,3]);
        se = strel('disk', 1);
        I_BW_e = imerode(I_BW_m, se);
        BWnobord = imclearborder(I_BW_e, 4);

        CC = bwconncomp(BWnobord, 8);
        stats = regionprops(CC, 'Centroid');
        centroids = cat(1, stats.Centroid);

        % Tally c-fos+ cells inside each ROI
        numROIs = length(roiData);
        regionIDs = zeros(numROIs,1);
        acronyms = cell(numROIs,1);
        regionNames = cell(numROIs,1);
        cellCounts = zeros(numROIs,1);
        areas = zeros(numROIs,1);

        for r = 1:numROIs
            % Try to parse region ID from ROI filename if you wish
            regionId = str2double(erase(roiData(r).fileName, '.roi'));
            regionIDs(r) = regionId;
            if ~isnan(regionId) && isKey(regionInfoMap, regionId)
                regionInfo = regionInfoMap(regionId);
                regionAcronym = regionInfo.acronym;
                regionName = regionInfo.name;
            else
                regionAcronym = 'Unknown';
                regionName = 'Unknown';
            end

            totalArea = 0;
            totalCount = 0;
            coords = roiData(r).numericCoords;

            if iscell(coords)
                % multiple polygons
                for c = 1:length(coords)
                    polyC = coords{c};
                    areaPoly = polyarea(polyC(:,1), polyC(:,2));
                    totalArea = totalArea + areaPoly;
                    inPoly = inpolygon(centroids(:,1), centroids(:,2), polyC(:,1), polyC(:,2));
                    totalCount = totalCount + sum(inPoly);
                end
            else
                % single polygon
                polyC = coords;
                totalArea = polyarea(polyC(:,1), polyC(:,2));
                inPoly = inpolygon(centroids(:,1), centroids(:,2), polyC(:,1), polyC(:,2));
                totalCount = sum(inPoly);
            end

            acronyms{r} = regionAcronym;
            regionNames{r} = regionName;
            cellCounts(r) = totalCount;
            areas(r) = totalArea;
        end

        %% 4e. Display result and save merged plot as PNG
        hFig = figure('Name', ['Merged TIF with ROIs: ' tifFiles(iFile).name], ...
                      'NumberTitle', 'off', 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        imshow(merged); hold on;
        for r = 1:numROIs
            coords = roiData(r).numericCoords;
            if isempty(coords), continue; end

            if iscell(coords)
                for c = 1:length(coords)
                    polyC = coords{c};
                    plot(polyC(:,1), polyC(:,2), 'r', 'LineWidth', 2);
                    centroid = mean(polyC,1);
                    regionId = str2double(erase(roiData(r).fileName, '.roi'));
                    if ~isnan(regionId) && isKey(regionInfoMap, regionId)
                        regionInfo = regionInfoMap(regionId);
                        text(centroid(1), centroid(2), regionInfo.acronym, ...
                            'Color', 'y', 'FontSize', 10, 'FontWeight', 'bold', ...
                            'HorizontalAlignment', 'center');
                    end
                end
            else
                plot(coords(:,1), coords(:,2), 'r', 'LineWidth', 2);
                centroid = mean(coords,1);
                regionId = str2double(erase(roiData(r).fileName, '.roi'));
                if ~isnan(regionId) && isKey(regionInfoMap, regionId)
                    regionInfo = regionInfoMap(regionId);
                    text(centroid(1), centroid(2), regionInfo.acronym, ...
                        'Color', 'y', 'FontSize', 10, 'FontWeight', 'bold', ...
                        'HorizontalAlignment', 'center');
                end
            end
        end
        hold off;
        
        % Save the figure as a PNG file with prefix "Merged_ROI_"
        [~, tifName, ~] = fileparts(tifFiles(iFile).name);
        pngName = fullfile(imageFolder, ['Merged_ROI_' tifName '.png']);
        saveas(hFig, pngName);
        fprintf('Merged plot saved: %s\n', pngName);
        close;

        %% 4f. Save summary
        summaryTable = table(regionIDs, acronyms, regionNames, cellCounts, areas, ...
            'VariableNames', {'ID','Acronym','Name','CellCount','Area_pixels'});
        excelName = fullfile(imageFolder, ['Summary of c-fos ' tifName '.xlsx']);
        writetable(summaryTable, excelName);
        fprintf('Summary saved: %s\n', excelName);

        % Optionally remove temp ROI folder contents:
        % rmdir(tempRoiDir, 's');

    catch ME
        fprintf(2, 'Error processing file %s: %s\n', tifFiles(iFile).name, ME.message);
    end
end

%% Local function: Recursively traverse JSON
function mapOut = traverseJson(node, mapIn)
    % Check that the node contains both acronym and name information
    if isfield(node, 'id') && isfield(node, 'acronym') && isfield(node, 'name')
        id = node.id;
        infoStruct.acronym = node.acronym;
        infoStruct.name = node.name;
        mapIn(id) = infoStruct;
    end
    if isfield(node, 'children') && ~isempty(node.children)
        for j = 1:length(node.children)
            mapIn = traverseJson(node.children(j), mapIn);
        end
    end
    mapOut = mapIn;
end

%% Local function: parse vfShapes
function polygons = parseVfShapes(vfShapes)
    if ischar(vfShapes)
        vf = sscanf(vfShapes, '%f');
    else
        vf = vfShapes;
    end
    polygons = {};
    currentPoly = [];
    i = 1;
    while i <= numel(vf)
        cmd = vf(i);
        i = i + 1;
        switch cmd
            case 0  % move-to
                if ~isempty(currentPoly)
                    polygons{end+1} = currentPoly;
                end
                currentPoly = [];
                if i+1 <= numel(vf)
                    x = vf(i); y = vf(i+1);
                    i = i + 2;
                    currentPoly = [x, y];
                end
            case 1  % line-to
                if i+1 <= numel(vf)
                    x = vf(i); y = vf(i+1);
                    i = i + 2;
                    currentPoly = [currentPoly; x, y];
                end
            case 4  % close-path
                if ~isempty(currentPoly)
                    if ~isequal(currentPoly(1,:), currentPoly(end,:))
                        currentPoly(end+1,:) = currentPoly(1,:);
                    end
                    polygons{end+1} = currentPoly;
                    currentPoly = [];
                end
            otherwise
                fprintf('Unrecognized command: %d\n', cmd);
        end
    end
    if ~isempty(currentPoly)
        polygons{end+1} = currentPoly;
    end
end
