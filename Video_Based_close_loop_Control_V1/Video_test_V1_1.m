%% Using videoinput with getsnapshot

% List available adaptors (e.g., 'winvideo')
hwInfo = imaqhwinfo;
disp(hwInfo.InstalledAdaptors);

% Create a video input object (adjust adaptor and device ID as needed)
vid = videoinput('winvideo', 1);
preview(vid);
disp('Live preview started. Press any key in the Command Window to capture a frame for ROI selection.');
pause;  % Wait for user input

% Close the preview
closepreview(vid);

% Capture a single frame using getsnapshot
frame = getsnapshot(vid);

% Display the captured frame
figure;
imshow(frame);
title('Captured Frame - Draw polygon to define ROI');

% Let user draw a polygon on the captured frame to define the ROI
hPoly = drawpolygon('LineWidth',2);
roi = hPoly.Position;  % This returns the polygon vertices
disp('Defined ROI coordinates:');
disp(roi);

% Optionally, save the captured frame to a file
imwrite(frame, 'captured_frame.png');
disp('Captured frame saved as "captured_frame.png".');

% Clean up the video input object
delete(vid);
clear vid;
