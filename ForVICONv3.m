clear
clc
close all

fps = 200; % 200 Hz = 1s
cd 'C:\Users\14000\Downloads\Year 4 FYP\Recruitment\Data Collection\P01_\Session 1\Snatch' % folder name was redacted
%%
viconfilename = 'P01_01_40_Snatch.csv';

vicondata = readmatrix(viconfilename, 'Range', 6); % Get whole Excel data from row 6
% Left1 marker
leftbar_x = vicondata(:,3)/1000; % Extract data and /1000 to convert to metre
leftbar_y = vicondata(:,4)/1000;
leftbar_z = vicondata(:,5)/1000; % Get vertical direction
% Right1 marker
rightbar_x = vicondata(:,15)/1000; % Extract data and convert to metre
rightbar_y = vicondata(:,16)/1000;
rightbar_z = vicondata(:,17)/1000; % Get vertical direction

leftbar_x(isnan(leftbar_x)) = 0; % Let NaN be 0
leftbar_y(isnan(leftbar_y)) = 0; % Let NaN be 0
leftbar_z(isnan(leftbar_z)) = 0; % Let NaN be 0
rightbar_x(isnan(rightbar_x)) = 0; % Let NaN be 0
rightbar_y(isnan(rightbar_y)) = 0; % Let NaN be 0
rightbar_z(isnan(rightbar_z)) = 0; % Let NaN be 0
% barpath_y = vicondata(:,4)/1000; % Get anterior posterior direction

% 15Hz low-pass Butterworth-Filter (4th order)
[d,c] = butter(4,15/(fps/2),'low');
leftbar_x_filter = filtfilt(d,c,leftbar_x); % filter
leftbar_y_filter = filtfilt(d,c,leftbar_y); % filter
leftbar_z_filter = filtfilt(d,c,leftbar_z); % filter
rightbar_x_filter = filtfilt(d,c,rightbar_x); % filter
rightbar_y_filter = filtfilt(d,c,rightbar_y); % filter
rightbar_z_filter = filtfilt(d,c,rightbar_z); % filter

% Create rep table
if str2double(viconfilename(8:9)) == 20 || str2double(viconfilename(8:9)) == 40 || str2double(viconfilename(8:9)) == 60
    repetition_number = array2table([1;2;3],'VariableNames',{'Repetition'});
elseif str2double(viconfilename(8:9)) == 70 || str2double(viconfilename(8:9)) == 80
   repetition_number = array2table([1;2],'VariableNames',{'Repetition'});
else
    repetition_number = array2table([1],'VariableNames',{'Repetition'});
end

% Getting time-series data for Force
starttime = 0; increment = 1/fps; timeend = length(leftbar_z_filter);
totaltime = starttime + increment*(0:timeend-1);

%% Left bar
% Detect start phase based on velocity graph
left_velocity_z = [];
for i=1:length(leftbar_z_filter)-1
    left_velocity_z(end+1) = (leftbar_z_filter(i+1)-leftbar_z_filter(i))/(1/fps);
end
% plot(leftbar_z_filter); hold on; plot(velocity_z); hold off;

leftvelocitystart = [];
for i=50:length(left_velocity_z)-50
    if left_velocity_z(i)>=0 && all(left_velocity_z(i:i+50) >= 0)
        leftvelocitystart(end+1) = i;
    end
end

% Remove excessive detection
for j=1:10
    try
        for i=1:length(leftvelocitystart)
            if i==length(leftvelocitystart)
                break
            end
            if ((leftvelocitystart(i+1) - leftvelocitystart(i)) < 200)
                leftvelocitystart(i+1) = [];
            end
        end
    catch
        continue
    end
end

% If need to manually remove start detection

% By index
leftindextoremove = [2]; % Check graph and change value inside []
leftvelocitystart(leftindextoremove) = [];

% Plot graph to see
figure('Name','Velocity Curve')
plot(left_velocity_z); hold on; plot(leftvelocitystart, left_velocity_z(leftvelocitystart), 'go'); 
legend('Velocity', 'Start'); hold off;

% figure('Name','Displacement Curve')
% plot(leftbar_z_filter); hold on; plot(leftvelocitystart, leftbar_z_filter(leftvelocitystart), 'go'); 
% legend('Displacement', 'Start'); hold off;
%% Find peak height
% Find peak displacement
[leftpeakheight, idx_leftpeakheight] = findpeaks(leftbar_z_filter, 'MinPeakHeight', 1, 'MinPeakDistance', 1000);

% Showing start to end
figure('Name','Displacement Curve');
plot(leftbar_z_filter); hold on; plot(leftvelocitystart, leftbar_z_filter(leftvelocitystart), 'go')
plot(idx_leftpeakheight, leftbar_z_filter(idx_leftpeakheight), 'ro'); 
legend('Displacement', 'Start', 'End'); hold off;
%%
% Plot each relevant left displacement curve
figure('Name','Displacement Curve - Each repetition');
for i=1:length(leftvelocitystart)
    plot(leftbar_z_filter(leftvelocitystart(i):idx_leftpeakheight(i)));
    hold on;
    left_trim_z{i} = (leftbar_z_filter(leftvelocitystart(i):idx_leftpeakheight(i)));
    if i==length(leftvelocitystart)
        hold off;
    end
end

% Plot bar path - left
figure('Name','Left Bar path - Each repetition');
for i=1:length(leftvelocitystart)
    plot(leftbar_x_filter(leftvelocitystart(i):idx_leftpeakheight(i)), leftbar_z_filter(leftvelocitystart(i):idx_leftpeakheight(i)))
    hold on;
%     left_trim_z{i} = (leftbar_z_filter(leftstart(i):idx_leftpeakheight(i)));
    if i==length(leftvelocitystart)
        hold off;
    end
end

% For each relevant cell, compute left velocity
for i=1:length(left_trim_z)
    left_velocity_z = [];
    for j=2:length(left_trim_z{i})
        left_velocity_z(end+1) = ((left_trim_z{i}(j))-(left_trim_z{i}(j-1)))/(1/fps);
    end
    rep_left_velocity_z{i} = left_velocity_z;
end

% Calculate left mean velocity
for i=1:length(rep_left_velocity_z)
    mean_left_velocity{i} = mean(rep_left_velocity_z{i});
end
mean_left_velocity = reshape(mean_left_velocity,[size(mean_left_velocity,2) 1]);
leftmeanvelocity = cell2table(mean_left_velocity, ...
"VariableNames",["LeftMeanVelocity"]);
leftmeanvelocity = [repetition_number leftmeanvelocity];
% Calculate left peak velocity
for i=1:length(rep_left_velocity_z)
    peak_left_velocity{i} = findpeaks(rep_left_velocity_z{i}, 'MinPeakHeight', 0.5, 'MinPeakDistance', 100);
end
peak_left_velocity = reshape(peak_left_velocity,[size(peak_left_velocity,2) 1]);
leftpeakvelocity = cell2table(peak_left_velocity, ...
"VariableNames",["LeftPeakVelocity"]);
leftpeakvelocity = [repetition_number leftpeakvelocity];
leftvelocity = join(leftmeanvelocity, leftpeakvelocity);

%% Right bar marker -----------------------------------------------------

% Detect start phase based on velocity graph
right_velocity_z = [];
for i=1:length(rightbar_z_filter)-1
    right_velocity_z(end+1) = (rightbar_z_filter(i+1)-rightbar_z_filter(i))/(1/fps);
end
% plot(leftbar_z_filter); hold on; plot(velocity_z); hold off;

rightvelocitystart = [];
for i=50:length(right_velocity_z)-50
    if right_velocity_z(i)>=0 && all(right_velocity_z(i:i+50) >= 0)
        rightvelocitystart(end+1) = i;
    end
end

% Remove excessive detection
for j=1:10
    try
        for i=1:length(rightvelocitystart)
            if i==length(rightvelocitystart)
                break
            end
            if ((rightvelocitystart(i+1) - rightvelocitystart(i)) < 200)
                rightvelocitystart(i+1) = [];
            end
        end
    catch
        continue
    end
end

% If need to manually remove start detection

% By index
rightindextoremove = [2]; % Check graph and change value inside []
rightvelocitystart(rightindextoremove) = [];

% Plot graph to see
% figure('Name','Velocity Curve')
% plot(right_velocity_z); hold on; plot(rightvelocitystart, right_velocity_z(rightvelocitystart), 'go'); 
% legend('Velocity', 'Start'); hold off;

% figure('Name','Displacement Curve')
% plot(rightbar_z_filter); hold on; plot(rightvelocitystart, rightbar_z_filter(rightvelocitystart), 'go'); 
% legend('Displacement', 'Start'); hold off;
%%
% Find peak displacement
[rightpeakheight, idx_rightpeakheight] = findpeaks(rightbar_z_filter, 'MinPeakHeight', 1, 'MinPeakDistance', 1000, 'NPeaks', 1);

% Showing start to end
figure('Name','Displacement Curve');
plot(rightbar_z_filter); hold on; plot(rightvelocitystart, rightbar_z_filter(rightvelocitystart), 'go')
plot(idx_rightpeakheight, rightbar_z_filter(idx_rightpeakheight), 'ro'); 
legend('Displacement', 'Start', 'End'); hold off;
%%
% Plot each relevant displacement curve
figure('Name','Displacement Curve - Each repetition');
for i=1:length(rightvelocitystart)
    plot(rightbar_z_filter(rightvelocitystart(i):idx_rightpeakheight(i)));
    hold on;
    right_trim_z{i} = (rightbar_z_filter(rightvelocitystart(i):idx_rightpeakheight(i)));
    if i==length(rightvelocitystart)
        hold off;
    end
end

% Plot bar path - right
figure('Name','Right Bar path - Each repetition');
for i=1:length(rightvelocitystart)
    plot(rightbar_x_filter(rightvelocitystart(i):idx_rightpeakheight(i)), rightbar_z_filter(rightvelocitystart(i):idx_rightpeakheight(i)))
    hold on;
%     left_trim_z{i} = (leftbar_z_filter(leftstart(i):idx_leftpeakheight(i)));
    if i==length(rightvelocitystart)
        hold off;
    end
end

% For each relevant cell, compute right velocity
for i=1:length(right_trim_z)
    right_velocity_z = [];
    for j=2:length(right_trim_z{i})
        right_velocity_z(end+1) = ((right_trim_z{i}(j))-(right_trim_z{i}(j-1)))/(1/fps);
    end
    rep_right_velocity_z{i} = right_velocity_z;
end

% Calculate mean right velocity
for i=1:length(rep_right_velocity_z)
    mean_right_velocity{i} = mean(rep_right_velocity_z{i});
end
mean_right_velocity = reshape(mean_right_velocity,[size(mean_right_velocity,2) 1]);
rightmeanvelocity = cell2table(mean_right_velocity, ...
"VariableNames",["RightMeanVelocity"]);
rightmeanvelocity = [repetition_number rightmeanvelocity];

% Calculate peak right velocity
for i=1:length(rep_right_velocity_z)
    peak_right_velocity{i} = findpeaks(rep_right_velocity_z{i}, 'MinPeakHeight', 0.7, 'MinPeakDistance', 100, 'NPeaks', 1);
end
peak_right_velocity = reshape(peak_right_velocity,[size(peak_right_velocity,2) 1]);
rightpeakvelocity = cell2table(peak_right_velocity, ...
"VariableNames",["RightPeakVelocity"]);
rightpeakvelocity = [repetition_number rightpeakvelocity];
rightvelocity = join(rightmeanvelocity, rightpeakvelocity);
%% Displacement curve of both bar ends
% Combine left and right bar path
figure('Name','Bar path - Each repetition');
barpath = tiledlayout(1,2);
nexttile;
for i=1:length(leftvelocitystart)
    txt = num2str(i);
    plot(leftbar_x_filter(leftvelocitystart(i):idx_leftpeakheight(i)), leftbar_z_filter(leftvelocitystart(i):idx_leftpeakheight(i)), 'DisplayName', txt)
    hold on; 
%     left_trim_z{i} = (leftbar_z_filter(leftstart(i):idx_leftpeakheight(i)));
    if i==length(leftvelocitystart)
        hold off;
        title('Left bar');
    end
end
legend show
nexttile;
for i=1:length(rightvelocitystart)
    txt = num2str(i);
    plot(rightbar_x_filter(rightvelocitystart(i):idx_rightpeakheight(i)), rightbar_z_filter(rightvelocitystart(i):idx_rightpeakheight(i)),  'DisplayName', txt)
    hold on;
%     left_trim_z{i} = (leftbar_z_filter(leftstart(i):idx_leftpeakheight(i)));
    if i==length(rightvelocitystart)
        hold off;
        title('Right bar');
    end
end
legend show
xlabel(barpath, 'Horizontal Displacement (m)')
ylabel(barpath, 'Vertical Displacement (m)')
%% Plot left and right, x and z axis
% Change the xlim and ylim where necessary
figure('Name', 'X-axis')
displacementfigure = tiledlayout(2,1);
nexttile;
plot(leftbar_x_filter); hold on; plot(rightbar_x_filter); hold off;
legend('Left', 'Right'); title('X - axis');
nexttile
plot(leftbar_z_filter); hold on; plot(rightbar_z_filter); hold off;
legend('Left', 'Right'); title('Z - axis');
xlabel(displacementfigure, 'No. of frames')
ylabel(displacementfigure, 'Displacement (m)')
%% Plot just 1
% figure('Name','Bar path - Each repetition');
% barpath = tiledlayout(1,2);
% nexttile;
% plot(leftbar_x_filter(leftstart(1):idx_leftpeakheight(1)), leftbar_z_filter(leftstart(1):idx_leftpeakheight(1)))
% title('Left bar'); xlim([-0.2 0.2]);
% nexttile;
% plot(rightbar_x_filter(rightstart(1):idx_rightpeakheight(1)), rightbar_z_filter(rightstart(1):idx_rightpeakheight(1)))
% title('Right bar'); xlim([-0.2 0.2]);
% xlabel(barpath, 'Horizontal Displacement (m)')
% ylabel(barpath, 'Vertical Displacement (m)')
%% Combine left and right velocity data
% Combine velocity data into 1 table
combinevelocitydata = join(leftvelocity, rightvelocity);

% Create %1RM table
% Get number of rows, check the %1RM based on filename, create an array of the same number, transpose and change to table
onerm = array2table(transpose(repelem((str2double(viconfilename(8:9))).',size(combinevelocitydata,1))), 'VariableName', {'%1RM'});
finaldata = [onerm combinevelocitydata];
disp(finaldata);
% Save to Excel

%% 
% writetable(finaldata, 'VelocityData.xlsx', 'WriteMode', 'Append')
