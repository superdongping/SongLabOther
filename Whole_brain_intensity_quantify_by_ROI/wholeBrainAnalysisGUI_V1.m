function wholeBrainAnalysisGUI_V1()
    % wholeBrainAnalysisGUI - Analyze serial TIFF images with 3 channels.
    % This modified version:
    %   1. Displays each merged TIFF image in full-screen.
    %   2. Allows the user to repeatedly draw ROIs.
    %   3. Prompts for the nucleus name before drawing each ROI.
    %   4. Computes ROI properties and saves cropped ROI images.
    %   5. Quantifies the area (in pixels) of green (Channel 2) and red (Channel 3)
    %      signals that are above a threshold determined by the display range.
    %   6. Optionally displays a demo of the thresholded images for pilot testing.
    %   7. Saves raw data (including nucleus name) to an Excel file and normalizes
    %      the signal areas across nuclei.
    clc;
    
    % Ask user if they want to run in demo mode.
    demoChoice = questdlg('Would you like to run in demo mode (show thresholded images)?', ...
        'Demo Mode', 'Yes', 'No', 'Yes');
    demoMode = strcmp(demoChoice, 'Yes');
    
    % Select folder containing TIFF files
    folderPath = uigetdir(pwd, 'Select folder with TIFF images');
    if folderPath == 0
        disp('No folder selected. Exiting.');
        return;
    end
    
    tifFiles = dir(fullfile(folderPath, '*.tif'));
    if isempty(tifFiles)
        error('No TIFF files found in the selected folder.');
    end
    
    % Initialize results table to store ROI data from all images.
    % Columns: FileName, ROI_Number, Nucleus, ROI_Area, Green_Area, Red_Area
    results = table('Size',[0 6],...
        'VariableTypes',{'string','double','string','double','double','double'},...
        'VariableNames',{'FileName','ROI_Number','Nucleus','ROI_Area','Green_Area','Red_Area'});
    
    % Process each TIFF file
    for i = 1:length(tifFiles)
        tifFile = fullfile(folderPath, tifFiles(i).name);
        info = imfinfo(tifFile);
        numPages = numel(info);
        if numPages < 3
            warning('File %s does not have at least 3 channels. Skipping.', tifFiles(i).name);
            continue;
        end
        
        % Read and store the image data (converted to double)
        img = zeros(info(1).Height, info(1).Width, numPages, 'double');
        for k = 1:numPages
            img(:,:,k) = double(imread(tifFile, k, 'Info', info));
        end
        
        % Process channels using specified intensity ranges
        % Channel 1 (Blue, DAPI): display range [0, 5000]
        blue = img(:,:,1);
        blue = min(max(blue, 0), 5000);
        blue_adjusted = blue / 5000;
        
        % Channel 2 (Green): display range [1000, 2000]
        green = img(:,:,2);
        green = min(max(green, 1000), 2000);
        green_adjusted = (green - 1000) / 1000;
        
        % Channel 3 (Red): display range [500, 2000]
        red = img(:,:,3);
        red = min(max(red, 500), 2000);
        red_adjusted = (red - 500) / 1500;
        
        % Create merged RGB image (order: [Red, Green, Blue])
        merged = cat(3, red_adjusted, green_adjusted, blue_adjusted);
        
        % Display merged image in full-screen
        hFig = figure('Name', sprintf('Analyzing: %s', tifFiles(i).name),...
                      'NumberTitle', 'off', 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        imshow(merged);
        title(sprintf('File: %s', tifFiles(i).name), 'Interpreter', 'none');
        drawnow;
        
        roiCount = 0;
        while true
            % Ask the user if they want to draw another ROI
            choice = questdlg('Do you want to draw another ROI?', ...
                'ROI Selection', 'Yes','No','Yes');
            if strcmp(choice, 'No')
                break;
            end
            
            % Prompt for nucleus name BEFORE drawing ROI
            prompt = {'Enter nucleus name for this ROI:'};
            dlg_title = 'Nucleus Name';
            answer = inputdlg(prompt, dlg_title, 1);
            if isempty(answer)
                nucName = "Unnamed";
            else
                nucName = string(answer{1});
            end
            
            % Let the user draw a freehand ROI
            hROI = drawfreehand('Color','g');
            roiMask = createMask(hROI);
            
            % Optionally display the ROI boundary on the image
            hold on;
            plot(hROI.Position(:,1), hROI.Position(:,2), 'LineWidth',2);
            hold off;
            
            roiCount = roiCount + 1;
            % Calculate ROI properties: area (in pixels)
            areaVal = sum(roiMask(:));
            
            % Quantify green signal area using a threshold method.
            % The threshold here (0.5) is chosen for the normalized green image.
            green_bw = imbinarize(green_adjusted, 0.5);
            green_areaROI = sum(green_bw(roiMask));
            
            % Quantify red signal area using a threshold method.
            red_bw = imbinarize(red_adjusted, 0.5);
            red_areaROI = sum(red_bw(roiMask));
            
            % If in demo mode, display the thresholded images for inspection.
            if demoMode
                hDemo = figure('Name', sprintf('Threshold Demo: ROI %d', roiCount));
                subplot(1,2,1);
                imshow(green_bw);
                hold on;
                % Highlight the ROI boundary in red for clarity.
                visboundaries(roiMask, 'Color', 'r', 'LineWidth', 1);
                title('Thresholded Green Channel');
                
                subplot(1,2,2);
                imshow(red_bw);
                hold on;
                visboundaries(roiMask, 'Color', 'r', 'LineWidth', 1);
                title('Thresholded Red Channel');
                % Wait for the user to close the demo window or press a key.
                uiwait(hDemo);
            end
            
            % Save the ROI image as a PNG file using the ROI bounding box
            stats = regionprops(roiMask, 'BoundingBox');
            if ~isempty(stats)
                bbox = round(stats(1).BoundingBox);
                % Ensure bounding box is within image bounds
                x1 = max(1, bbox(1));
                y1 = max(1, bbox(2));
                x2 = min(size(merged,2), x1 + bbox(3));
                y2 = min(size(merged,1), y1 + bbox(4));
                roiImage = merged(y1:y2, x1:x2, :);
                [~, baseName, ~] = fileparts(tifFiles(i).name);
                pngName = sprintf('%s_ROI%d_%s.png', baseName, roiCount, nucName);
                imwrite(roiImage, fullfile(folderPath, pngName));
            end
            
            % Append ROI results to the table (including nucleus name)
            newRow = {tifFiles(i).name, roiCount, nucName, areaVal, green_areaROI, red_areaROI};
            results = [results; newRow]; %#ok<AGROW>
            
            % Delete ROI object for next iteration
            delete(hROI);
        end
        
        % Close the full-screen figure for the current image
        close(hFig);
    end  % End for each TIFF file
    
    % Save raw data to Excel (including nucleus name)
    excelFile = fullfile(folderPath, 'Summary_output.xlsx');
    writetable(results, excelFile, 'Sheet', 'Raw_Data');
    
    % --- Normalized Data ---
    % Group results by nucleus name
    [groupIDs, groupNames] = findgroups(results.Nucleus);
    
    % For each nucleus group, sum the Green_Area and Red_Area.
    combinedGreen = splitapply(@sum, results.Green_Area, groupIDs);
    combinedRed   = splitapply(@sum, results.Red_Area, groupIDs);
    
    % Calculate total areas across all nuclei for each channel.
    totalGreenAll = sum(combinedGreen);
    totalRedAll   = sum(combinedRed);
    
    % Normalize each nucleus's channel separately.
    normalizedGreen = (combinedGreen / totalGreenAll) * 100;
    normalizedRed   = (combinedRed / totalRedAll) * 100;
    
    % Create a table for normalized data.
    normData = table(groupNames, combinedGreen, combinedRed, normalizedGreen, normalizedRed, ...
        'VariableNames',{'Nucleus','Combined_Green_Area','Combined_Red_Area','Normalized_Green','Normalized_Red'});
    
    % Save normalized data to the same Excel file in a new sheet.
    writetable(normData, excelFile, 'Sheet', 'Normalized_Data');
    
    disp(['Analysis complete. Summary saved to ' excelFile]);
end
