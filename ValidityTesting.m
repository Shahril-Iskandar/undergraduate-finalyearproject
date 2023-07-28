clear
clc
close all
cd 'C:\Users\14000\Downloads\Year 4 FYP\Recruitment\Data Collection'
%%
powerclean2 = readmatrix('Participant Info.xlsx', 'Sheet', 'D1 - Power Clean (combine)');
snatch2 = readmatrix('Participant Info.xlsx', 'Sheet', 'D1 - Snatch (Combine)');

powerclean_meandata = readmatrix('Dataset.xlsx', 'Sheet', 'PowerClean_Mean');
powerclean_peakdata = readmatrix('Dataset.xlsx', 'Sheet', 'PowerClean_Peak');
snatch_meandata = readmatrix('Dataset.xlsx', 'Sheet', 'Snatch_Mean');
snatch_peakdata = readmatrix('Dataset.xlsx', 'Sheet', 'Snatch_Peak');
intensity = [20 40 60 70 80 90];

% Remove rows with NaN (blank cells)
powerclean_meandata = powerclean_meandata(sum(isnan(powerclean_meandata),2)==0,:);
powerclean_peakdata = powerclean_peakdata(sum(isnan(powerclean_peakdata),2)==0,:);
snatch_meandata = snatch_meandata(sum(isnan(snatch_meandata),2)==0,:);
snatch_peakdata = snatch_peakdata(sum(isnan(snatch_peakdata),2)==0,:);

% Add average of the 2 devices as Column 6
powerclean_meandata = [powerclean_meandata (powerclean_meandata(:,4)+powerclean_meandata(:,5))/2];
powerclean_peakdata = [powerclean_peakdata (powerclean_peakdata(:,4)+powerclean_peakdata(:,5))/2];
snatch_meandata = [snatch_meandata (snatch_meandata(:,4)+snatch_meandata(:,5))/2];
snatch_peakdata = [snatch_peakdata (snatch_peakdata(:,4)+snatch_peakdata(:,5))/2];

% Add differences of the 2 devices as Column 7
powerclean_meandata = [powerclean_meandata (powerclean_meandata(:,4)-powerclean_meandata(:,5))];
powerclean_peakdata = [powerclean_peakdata (powerclean_peakdata(:,4)-powerclean_peakdata(:,5))];
snatch_meandata = [snatch_meandata (snatch_meandata(:,4)-snatch_meandata(:,5))];
snatch_peakdata = [snatch_peakdata (snatch_peakdata(:,4)-snatch_peakdata(:,5))];

% sortrows(powerclean2,1) % Sort rows in column 3

% Power Clean Mean Bland Altman
for i=1:length(intensity)
    data = powerclean_meandata(powerclean_meandata(:,2) == intensity(i),:); % Get values of specific intensity according to loop
    meandifferences = mean(data(:,7)); % Calculate the mean of the differences
    lowersd = meandifferences-1.96*std(data(:,7)); % Calculate lower LoA
    uppersd = meandifferences+1.96*std(data(:,7)); % Calculate upper LoA
    sprintf('The mean for %d%% is %f m/s', intensity(i), meandifferences)
    sprintf('Limits of agreement for %d%% is %f m/s to %f m/s', intensity(i), lowersd, uppersd)
    slope = data(:,4)\data(:,5);
    sprintf('Slope / Regression coefficient: %f', slope);
%     figure('Name',sprintf('%d%% 1RM Bland-Altman Plot', intensity(i)));
    subplot(2,3,i);
    scatter(data(:,6), data(:,7), 'filled', 'MarkerFaceColor', [0 0 1]);  % Plot mean velocity Bland Altman 
%     hold on;
    yline(uppersd, '--r');
    % textscatter(1.21, left_upper_40+0.01, "+1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
    % textscatter(1.21, left_upper_40-0.01, string(num2str(left_upper_40)), 'ColorData', [1 0 0])% Change first value accordingly
    yline(meandifferences, 'b');
    % textscatter(1.21, left_average_difference_40+0.01, "Mean", 'ColorData', [0 0 1]) % Change first value accordingly
    % textscatter(1.21, left_average_difference_40-0.01, string(num2str(left_average_difference_40)), 'ColorData', [0 0 1]) % Change first value accordingly
    yline(lowersd, '--r');
    % textscatter(1.21, left_lower_40+0.01, "-1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
    % textscatter(1.21, left_lower_40-0.01, string(num2str(left_lower_40)), 'ColorData', [1 0 0]) % Change first value accordingly
    % axis equal
    lsline();
    title(sprintf('%d%% 1RM Bland-Altman Plot', intensity(i))); xlabel('Mean of FLEX and VICON'); ylabel('Difference between FLEX and VICON');
end

%%
% Read participant info excel, sheet D1 power clean
d1powerclean = readmatrix('Participant Info.xlsx', 'Sheet', 'D1 - Power Clean'); % D1 - Power Clean
d2powerclean = readmatrix('Participant Info.xlsx', 'Sheet', 'D2 - Power Clean'); % D1 - Power Clean
d1snatch= readmatrix('Participant Info.xlsx', 'Sheet', 'D1 - Snatch'); % D1 - Power Clean
d2snatch = readmatrix('Participant Info.xlsx', 'Sheet', 'D2 - Snatch'); % D1 - Power Clean

% d1powerclean_flex_meanvelocity_left = d1powerclean(:,4); % D1 FLEX_Mean_Concentric_Left
% d1powerclean_flex_peakvelocity_left = d1powerclean(:,5); % D1 FLEX_Peak_Concentric_Left
% d1powerclean_flex_meanvelocity_left = d1powerclean(:,6); % D1 FLEX_Mean_Concentric_Right
% d1powerclean_flex_peakvelocity_left = d1powerclean(:,7); % D1 FLEX_Peak_Concentric_Right

% Combine all data
powerclean = cat(1,d1powerclean,d2powerclean);
snatch = cat(1,d1snatch,d2snatch);
alldata = cat(1,powerclean,snatch);

% Remove rows with NaN
alldata(any(isnan(alldata),2),:) = [];

% Move everything alldata_xx up.
% Separate left and right
% Analyse

% Column 4 = FLEX Mean Left
% Column 5 = FLEX Peak Left
% Column 6 = FLEX Mean Right
% Column 7 = FLEX Peak Right
% Column 8 = VICON Mean Left
% Column 9 = VICON Peak Left
% Column 10 = VICON Mean Right
% Column 11 = VICON Peak Right
% To compare Mean Left: Column 4 and 8
% To compare Peak Left: Column 5 and 9
% To compare Mean Right: Column 6 and 10
% To compare Peak Right: Column 7 and 11
%% Extract only 40% 1RM data
data40 = alldata(:,2) == 40; % Logical array to find if column 2 has 40
alldata_40 = alldata(data40,:);

% Get FLEX and VICON mean concentric (Left)
left_trim_40 = [alldata_40(:,4) alldata_40(:,8)];
left_average_40 = mean(left_trim_40,2);
left_difference_40 = left_trim_40(:,1) - left_trim_40(:,2);

% [H, pValue, W] = swtest(left_difference_40, 0.05)
% 
% if pValue <= 0.05
%     left_difference_40 = log(left_difference_40);
% end

left_average_difference_40 = mean(left_difference_40);
left_lower_40 = left_average_difference_40-1.96*std(left_difference_40);
left_upper_40 = left_average_difference_40+1.96*std(left_difference_40);
% see_40 = sum(left_trim_40(:,1) - left_trim_40(:,2)) / (size(left_trim_40,1) - 2); % Standard error
% cv_40 = std(left_trim_40(:,1))/mean(left_trim_40(:,1));
% cv_40vicon = std(left_trim_40(:,2))/mean(left_trim_40(:,2));

% Plot Bland-Altman graph (Left)
figure('Name', '40% 1RM Bland-Altman Plot (Left)')
tiledlayout(1,2);
nexttile
scatter(left_average_40, left_difference_40, 'filled', 'MarkerFaceColor', [0 0 1]); hold on;
yline(left_upper_40, '--r');
textscatter(1.21, left_upper_40+0.01, "+1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
textscatter(1.21, left_upper_40-0.01, string(num2str(left_upper_40)), 'ColorData', [1 0 0])% Change first value accordingly
yline(left_average_difference_40, 'b');
textscatter(1.21, left_average_difference_40+0.01, "Mean", 'ColorData', [0 0 1]) % Change first value accordingly
textscatter(1.21, left_average_difference_40-0.01, string(num2str(left_average_difference_40)), 'ColorData', [0 0 1]) % Change first value accordingly
yline(left_lower_40, '--r');
textscatter(1.21, left_lower_40+0.01, "-1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
textscatter(1.21, left_lower_40-0.01, string(num2str(left_lower_40)), 'ColorData', [1 0 0]) % Change first value accordingly
lsline();
% axis equal
title('40% 1RM Bland-Altman Plot (Left)'); xlabel('Mean of FLEX and VICON'); ylabel('Difference between FLEX and VICON'); hold off;

% Plot correlation graph left
nexttile
% figure('Name', '40% 1RM Correlation Plot (Left)')
scatter(left_trim_40(:,1), left_trim_40(:,2), 'filled');
lsline(); refline(1);
title('40% 1RM Correlation Plot (Left)'); xlabel('FLEX - Mean Concentric Velocity (m/s)'); ylabel('VICON - Mean Concentric Velocity (m/s)')

left_R_40 = corrcoef(left_trim_40);
sprintf('40%% Pearsons Correlation r (Left): %f', left_R_40(2))

% Get FLEX and VICON mean concentric (Right)
right_trim_40 = [alldata_40(:,6) alldata_40(:,10)];
right_average_40 = mean(right_trim_40,2);
right_difference_40 = right_trim_40(:,1) - right_trim_40(:,2);
[H, pValue, W] = swtest(right_difference_40, 0.05)
right_average_difference_40 = mean(right_difference_40);
right_lower_40 = right_average_difference_40-1.96*std(right_difference_40);
right_upper_40 = right_average_difference_40+1.96*std(right_difference_40);
see_40 = sum(right_trim_40(:,1) - right_trim_40(:,2)) / (size(right_trim_40,1) - 2); % Standard error

% Plot Bland-Altman graph (Right)
figure('Name', '40% 1RM Bland-Altman Plot (Right)')
tiledlayout(1,2);
nexttile
scatter(right_average_40, right_difference_40, 'filled', 'MarkerFaceColor', [0 0 1]); hold on;
yline(right_upper_40, '--r');
textscatter(1.21, right_upper_40+0.01, "+1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
textscatter(1.21, right_upper_40-0.01, string(num2str(right_upper_40)), 'ColorData', [1 0 0])% Change first value accordingly
yline(right_average_difference_40, 'b');
textscatter(1.21, right_average_difference_40+0.01, "Mean", 'ColorData', [0 0 1]) % Change first value accordingly
textscatter(1.21, right_average_difference_40-0.01, string(num2str(right_average_difference_40)), 'ColorData', [0 0 1]) % Change first value accordingly
yline(right_lower_40, '--r');
textscatter(1.21, right_lower_40+0.01, "-1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
textscatter(1.21, right_lower_40-0.01, string(num2str(right_lower_40)), 'ColorData', [1 0 0]) % Change first value accordingly
title('40% 1RM Bland-Altman Plot (Right)'); xlabel('Mean of FLEX and VICON'); ylabel('Difference between FLEX and VICON'); hold off;

% Plot correlation graph right
nexttile
% figure('Name', '40% 1RM Correlation Plot (Left)')
scatter(right_trim_40(:,1), right_trim_40(:,2), 'filled');
lsline();
title('40% 1RM Correlation Plot (Right)'); xlabel('FLEX - Mean Concentric Velocity (m/s)'); ylabel('VICON - Mean Concentric Velocity (m/s)')

right_R_40 = corrcoef(right_trim_40);
sprintf('40%% Pearsons Correlation r (Right): %f', right_R_40(2))
%% Extract only 60% 1RM data
data60 = alldata(:,2) == 60;
alldata_60 = alldata(data60,:);

% Get FLEX and VICON mean concentric (Left)
left_trim_60 = [alldata_60(:,4) alldata_60(:,8)];
left_average_60 = mean(left_trim_60,2);
left_difference_60 = left_trim_60(:,1) - left_trim_60(:,2);
[H, pValue, W] = swtest(left_difference_60, 0.05)
left_average_difference_60 = mean(left_difference_60);
left_lower_60 = left_average_difference_60-1.96*std(left_difference_60);
left_upper_60 = left_average_difference_60+1.96*std(left_difference_60);

% Plot Bland-Altman graph
figure('Name', '60% 1RM Bland-Altman Plot (Left)')
tiledlayout(1,2);
nexttile;
scatter(left_average_60, left_difference_60, 'filled', 'MarkerFaceColor', [0 0 1]); hold on;
yline(left_upper_60, '--r');
textscatter(1.21, left_upper_60+0.01, "+1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
textscatter(1.21, left_upper_60-0.01, string(num2str(left_upper_60)), 'ColorData', [1 0 0])% Change first value accordingly
yline(left_average_difference_60, 'b');
textscatter(1.21, left_average_difference_60+0.01, "Mean", 'ColorData', [0 0 1]) % Change first value accordingly
textscatter(1.21, left_average_difference_60-0.01, string(num2str(left_average_difference_60)), 'ColorData', [0 0 1]) % Change first value accordingly
yline(left_lower_60, '--r');
textscatter(1.21, left_lower_60+0.01, "-1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
textscatter(1.21, left_lower_60-0.01, string(num2str(left_lower_60)), 'ColorData', [1 0 0]) % Change first value accordingly
title('60% 1RM Bland-Altman Plot (Left)'); xlabel('Mean of FLEX and VICON'); ylabel('Difference between FLEX and VICON'); hold off;

% Plot correlation graph
nexttile;
% figure('Name', '60% 1RM Correlation Plot (Left)')
scatter(left_trim_60(:,1), left_trim_60(:,2), 'filled');
lsline();
title('60% 1RM Correlation Plot (Left)'); xlabel('FLEX - Mean Concentric Velocity (m/s)'); ylabel('VICON - Mean Concentric Velocity (m/s)')

left_R_60 = corrcoef(left_trim_60);
sprintf('60%% Pearsons Correlation r (Left): %f', left_R_60(2))

% Get FLEX and VICON mean concentric (Right)
right_trim_60 = [alldata_60(:,6) alldata_60(:,10)];
right_average_60 = mean(right_trim_60,2);
right_difference_60 = right_trim_60(:,1) - right_trim_60(:,2);
[H, pValue, W] = swtest(right_difference_60, 0.05)
right_average_difference_60 = mean(right_difference_60);
right_lower_60 = right_average_difference_60-1.96*std(right_difference_60);
right_upper_60 = right_average_difference_60+1.96*std(right_difference_60);

% Plot Bland-Altman graph
figure('Name', '60% 1RM Bland-Altman Plot (Right)')
tiledlayout(1,2);
nexttile;
scatter(right_average_60, right_difference_60, 'filled', 'MarkerFaceColor', [0 0 1]); hold on;
yline(right_upper_60, '--r');
textscatter(1.21, right_upper_60+0.01, "+1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
textscatter(1.21, right_upper_60-0.01, string(num2str(right_upper_60)), 'ColorData', [1 0 0])% Change first value accordingly
yline(right_average_difference_60, 'b');
textscatter(1.21, right_average_difference_60+0.01, "Mean", 'ColorData', [0 0 1]) % Change first value accordingly
textscatter(1.21, right_average_difference_60-0.01, string(num2str(right_average_difference_60)), 'ColorData', [0 0 1]) % Change first value accordingly
yline(right_lower_60, '--r');
textscatter(1.21, right_lower_60+0.01, "-1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
textscatter(1.21, right_lower_60-0.01, string(num2str(right_lower_60)), 'ColorData', [1 0 0]) % Change first value accordingly
title('60% 1RM Bland-Altman Plot (Right)'); xlabel('Mean of FLEX and VICON'); ylabel('Difference between FLEX and VICON'); hold off;

% Plot correlation graph
nexttile;
% figure('Name', '60% 1RM Correlation Plot (Right)')
scatter(right_trim_60(:,1), right_trim_60(:,2), 'filled');
lsline();
title('60% 1RM Correlation Plot (Right)'); xlabel('FLEX - Mean Concentric Velocity (m/s)'); ylabel('VICON - Mean Concentric Velocity (m/s)')

right_R_60 = corrcoef(right_trim_60);
sprintf('60%% Pearsons Correlation r (Right): %f', right_R_60(2))
%% Extract only 70% 1RM data
data70 = alldata(:,2) == 70;
alldata_70 = alldata(data70,:);

% Get FLEX and VICON mean concentric (Left)
left_trim_70 = [alldata_70(:,4) alldata_70(:,8)];
left_average_70 = mean(left_trim_70,2);
left_difference_70 = left_trim_70(:,1) - left_trim_70(:,2);
[H, pValue, W] = swtest(left_difference_70, 0.05)
left_average_difference_70 = mean(left_difference_70);
left_lower_70 = left_average_difference_70-1.96*std(left_difference_70);
left_upper_70 = left_average_difference_70+1.96*std(left_difference_70);

% Plot Bland-Altman graph
figure('Name', '70% 1RM Bland-Altman Plot (Left)')
tiledlayout(1,2);
nexttile;
scatter(left_average_70, left_difference_70, 'filled', 'MarkerFaceColor', [0 0 1]); hold on;
yline(left_upper_70, '--r');
textscatter(1.21, left_upper_70+0.01, "+1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
textscatter(1.21, left_upper_70-0.01, string(num2str(left_upper_70)), 'ColorData', [1 0 0])% Change first value accordingly
yline(left_average_difference_70, 'b');
textscatter(1.21, left_average_difference_70+0.01, "Mean", 'ColorData', [0 0 1]) % Change first value accordingly
textscatter(1.21, left_average_difference_70-0.01, string(num2str(left_average_difference_70)), 'ColorData', [0 0 1]) % Change first value accordingly
yline(left_lower_70, '--r');
textscatter(1.21, left_lower_70+0.01, "-1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
textscatter(1.21, left_lower_70-0.01, string(num2str(left_lower_70)), 'ColorData', [1 0 0]) % Change first value accordingly
title('70% 1RM Bland-Altman Plot (Left)'); xlabel('Mean of FLEX and VICON'); ylabel('Difference between FLEX and VICON'); hold off;

% Plot correlation graph
nexttile;
% figure('Name', '70% 1RM Correlation Plot (Left)')
scatter(left_trim_70(:,1), left_trim_70(:,2), 'filled');
lsline();
title('70% 1RM Correlation Plot (Left)'); xlabel('FLEX - Mean Concentric Velocity (m/s)'); ylabel('VICON - Mean Concentric Velocity (m/s)')

left_R_70 = corrcoef(left_trim_70);
sprintf('70%% Pearsons Correlation r (Left): %f', left_R_70(2))

% Get FLEX and VICON mean concentric (Right)
right_trim_70 = [alldata_70(:,6) alldata_70(:,10)];
right_average_70 = mean(right_trim_70,2);
right_difference_70 = right_trim_70(:,1) - right_trim_70(:,2);
[H, pValue, W] = swtest(right_difference_70, 0.05)
right_average_difference_70 = mean(right_difference_70);
right_lower_70 = right_average_difference_70-1.96*std(right_difference_70);
right_upper_70 = right_average_difference_70+1.96*std(right_difference_70);

% Plot Bland-Altman graph
figure('Name', '70% 1RM Bland-Altman Plot (Right)')
tiledlayout(1,2);
nexttile;
scatter(right_average_70, right_difference_70, 'filled', 'MarkerFaceColor', [0 0 1]); hold on;
yline(right_upper_70, '--r');
textscatter(1.21, right_upper_70+0.01, "+1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
textscatter(1.21, right_upper_70-0.01, string(num2str(right_upper_70)), 'ColorData', [1 0 0])% Change first value accordingly
yline(right_average_difference_70, 'b');
textscatter(1.21, right_average_difference_70+0.01, "Mean", 'ColorData', [0 0 1]) % Change first value accordingly
textscatter(1.21, right_average_difference_70-0.01, string(num2str(right_average_difference_70)), 'ColorData', [0 0 1]) % Change first value accordingly
yline(right_lower_70, '--r');
textscatter(1.21, right_lower_70+0.01, "-1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
textscatter(1.21, right_lower_70-0.01, string(num2str(right_lower_70)), 'ColorData', [1 0 0]) % Change first value accordingly
title('70% 1RM Bland-Altman Plot (Right)'); xlabel('Mean of FLEX and VICON'); ylabel('Difference between FLEX and VICON'); hold off;

% Plot correlation graph
nexttile;
% figure('Name', '70% 1RM Correlation Plot (Right)')
scatter(right_trim_70(:,1), right_trim_70(:,2), 'filled');
lsline();
title('70% 1RM Correlation Plot (Right)'); xlabel('FLEX - Mean Concentric Velocity (m/s)'); ylabel('VICON - Mean Concentric Velocity (m/s)')

right_R_70= corrcoef(right_trim_70);
sprintf('70%% Pearsons Correlation r (Right): %f', right_R_70(2))
%% Extract only 80% 1RM data
data80 = alldata(:,2) == 80;
alldata_80 = alldata(data80,:);

% Get FLEX and VICON mean concentric (Left)
left_trim_80 = [alldata_80(:,4) alldata_80(:,8)];
left_average_80 = mean(left_trim_80,2);
left_difference_80 = left_trim_80(:,1) - left_trim_80(:,2);
[H, pValue, W] = swtest(left_difference_80, 0.05)
left_average_difference_80 = mean(left_difference_80);
left_lower_80 = left_average_difference_80-1.96*std(left_difference_80);
left_upper_80 = left_average_difference_80+1.96*std(left_difference_80);

% Plot Bland-Altman graph
figure('Name', '80% 1RM Bland-Altman Plot (Left)')
tiledlayout(1,2);
nexttile;
scatter(left_average_80, left_difference_80, 'filled', 'MarkerFaceColor', [0 0 1]); hold on;
yline(left_upper_80, '--r');
textscatter(1.21, left_upper_80+0.01, "+1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
textscatter(1.21, left_upper_80-0.01, string(num2str(left_upper_80)), 'ColorData', [1 0 0])% Change first value accordingly
yline(left_average_difference_80, 'b');
textscatter(1.21, left_average_difference_80+0.01, "Mean", 'ColorData', [0 0 1]) % Change first value accordingly
textscatter(1.21, left_average_difference_80-0.01, string(num2str(left_average_difference_80)), 'ColorData', [0 0 1]) % Change first value accordingly
yline(left_lower_80, '--r');
textscatter(1.21, left_lower_80+0.01, "-1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
textscatter(1.21, left_lower_80-0.01, string(num2str(left_lower_80)), 'ColorData', [1 0 0]) % Change first value accordingly
title('80% 1RM Bland-Altman Plot (Left)'); xlabel('Mean of FLEX and VICON'); ylabel('Difference between FLEX and VICON'); hold off;

% Plot correlation graph
nexttile;
% figure('Name', '80% 1RM Correlation Plot (Left)')
scatter(left_trim_80(:,1), left_trim_80(:,2), 'filled');
lsline();
title('80% 1RM Correlation Plot (Left)'); xlabel('FLEX - Mean Concentric Velocity (m/s)'); ylabel('VICON - Mean Concentric Velocity (m/s)')

left_R_80 = corrcoef(left_trim_80);
sprintf('80%% Pearsons Correlation r (Left): %f', left_R_80(2))

% Get FLEX and VICON mean concentric (Right)
right_trim_80 = [alldata_80(:,6) alldata_80(:,10)];
right_average_80 = mean(right_trim_80,2);
right_difference_80 = right_trim_80(:,1) - right_trim_80(:,2);
[H, pValue, W] = swtest(right_difference_80, 0.05)
right_average_difference_80 = mean(right_difference_80);
right_lower_80 = right_average_difference_80-1.96*std(right_difference_80);
right_upper_80 = right_average_difference_80+1.96*std(right_difference_80);

% Plot Bland-Altman graph
figure('Name', '80% 1RM Bland-Altman Plot (Right)')
tiledlayout(1,2);
nexttile;
scatter(right_average_80, right_difference_80, 'filled', 'MarkerFaceColor', [0 0 1]); hold on;
yline(right_upper_80, '--r');
textscatter(1.21, right_upper_80+0.01, "+1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
textscatter(1.21, right_upper_80-0.01, string(num2str(right_upper_80)), 'ColorData', [1 0 0])% Change first value accordingly
yline(right_average_difference_80, 'b');
textscatter(1.21, right_average_difference_80+0.01, "Mean", 'ColorData', [0 0 1]) % Change first value accordingly
textscatter(1.21, right_average_difference_80-0.01, string(num2str(right_average_difference_80)), 'ColorData', [0 0 1]) % Change first value accordingly
yline(right_lower_80, '--r');
textscatter(1.21, right_lower_80+0.01, "-1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
textscatter(1.21, right_lower_80-0.01, string(num2str(right_lower_80)), 'ColorData', [1 0 0]) % Change first value accordingly
title('80% 1RM Bland-Altman Plot (Right)'); xlabel('Mean of FLEX and VICON'); ylabel('Difference between FLEX and VICON'); hold off;

% Plot correlation graph
nexttile;
% figure('Name', '80% 1RM Correlation Plot (Left)')
scatter(right_trim_80(:,1), right_trim_80(:,2), 'filled');
lsline();
title('80% 1RM Correlation Plot (Right)'); xlabel('FLEX - Mean Concentric Velocity (m/s)'); ylabel('VICON - Mean Concentric Velocity (m/s)')

right_R_80 = corrcoef(right_trim_80);
sprintf('80%% Pearsons Correlation r (Right): %f', right_R_80(2))
%% Extract only 90% 1RM data
data90 = alldata(:,2) == 90;
alldata_90 = alldata(data90,:);

% Get FLEX and VICON mean concentric (Left)
left_trim_90 = [alldata_90(:,4) alldata_90(:,8)];
left_average_90 = mean(left_trim_90,2);
left_difference_90 = left_trim_90(:,1) - left_trim_90(:,2);
[H, pValue, W] = swtest(left_difference_80, 0.05)
left_average_difference_90 = mean(left_difference_90);
left_lower_90 = left_average_difference_90-1.96*std(left_difference_90);
left_upper_90 = left_average_difference_90+1.96*std(left_difference_90);

% Plot Bland-Altman graph
figure('Name', '90% 1RM Bland-Altman Plot (Left)')
tiledlayout(1,2);
nexttile;
scatter(left_average_90, left_difference_90, 'filled', 'MarkerFaceColor', [0 0 1]); hold on;
yline(left_upper_90, '--r');
textscatter(1.21, left_upper_90+0.01, "+1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
textscatter(1.21, left_upper_90-0.01, string(num2str(left_upper_90)), 'ColorData', [1 0 0])% Change first value accordingly
yline(left_average_difference_90, 'b');
textscatter(1.21, left_average_difference_90+0.01, "Mean", 'ColorData', [0 0 1]) % Change first value accordingly
textscatter(1.21, left_average_difference_90-0.01, string(num2str(left_average_difference_90)), 'ColorData', [0 0 1]) % Change first value accordingly
yline(left_lower_90, '--r');
textscatter(1.21, left_lower_90+0.01, "-1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
textscatter(1.21, left_lower_90-0.01, string(num2str(left_lower_90)), 'ColorData', [1 0 0]) % Change first value accordingly
title('90% 1RM Bland-Altman Plot (Left)'); xlabel('Mean of FLEX and VICON'); ylabel('Difference between FLEX and VICON'); hold off;

% Plot correlation graph
nexttile;
% figure('Name', '90% 1RM Correlation Plot (Left)')
scatter(left_trim_90(:,1), left_trim_90(:,2), 'filled');
lsline();
title('90% 1RM Correlation Plot (Left)'); xlabel('FLEX - Mean Concentric Velocity (m/s)'); ylabel('VICON - Mean Concentric Velocity (m/s)')

left_R_90 = corrcoef(left_trim_90);
sprintf('90%% Pearsons Correlation r (Left): %f', left_R_90(2))

% Get FLEX and VICON mean concentric (Right)
right_trim_90 = [alldata_90(:,6) alldata_90(:,10)];
right_average_90 = mean(right_trim_90,2);
right_difference_90 = right_trim_90(:,1) - right_trim_90(:,2);
[H, pValue, W] = swtest(right_difference_90, 0.05);
right_average_difference_90 = mean(right_difference_90);
right_lower_90 = right_average_difference_90-1.96*std(right_difference_90);
right_upper_90 = right_average_difference_90+1.96*std(right_difference_90);

% Plot Bland-Altman graph
figure('Name', '90% 1RM Bland-Altman Plot (Right)')
tiledlayout(1,2);
nexttile;
scatter(right_average_90, right_difference_90, 'filled', 'MarkerFaceColor', [0 0 1]); hold on;
yline(right_upper_90, '--r');
textscatter(1.21, right_upper_90+0.01, "+1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
textscatter(1.21, right_upper_90-0.01, string(num2str(right_upper_90)), 'ColorData', [1 0 0])% Change first value accordingly
yline(right_average_difference_90, 'b');
textscatter(1.21, right_average_difference_90+0.01, "Mean", 'ColorData', [0 0 1]) % Change first value accordingly
textscatter(1.21, right_average_difference_90-0.01, string(num2str(right_average_difference_90)), 'ColorData', [0 0 1]) % Change first value accordingly
yline(right_lower_90, '--r');
textscatter(1.21, right_lower_90+0.01, "-1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
textscatter(1.21, right_lower_90-0.01, string(num2str(right_lower_90)), 'ColorData', [1 0 0]) % Change first value accordingly
title('90% 1RM Bland-Altman Plot (Right)'); xlabel('Mean of FLEX and VICON'); ylabel('Difference between FLEX and VICON'); hold off;

% Plot correlation graph
nexttile;
% figure('Name', '90% 1RM Correlation Plot (Right)')
scatter(right_trim_90(:,1), right_trim_90(:,2), 'filled');
lsline();
title('90% 1RM Correlation Plot (Right)'); xlabel('FLEX - Mean Concentric Velocity (m/s)'); ylabel('VICON - Mean Concentric Velocity (m/s)')

right_R_90 = corrcoef(right_trim_90);
sprintf('90%% Pearsons Correlation r (Right): %f', right_R_90(2))
%% Combine all left data
% Get FLEX and VICON mean concentric
all_left = [alldata(:,4) alldata(:,8)];
left_average = mean(all_left,2);
left_difference = all_left(:,1) - all_left(:,2);
left_average_difference = mean(left_difference);
left_lower = left_average_difference-1.96*std(left_difference);
left_upper = left_average_difference+1.96*std(left_difference);

% Plot Bland-Altman graph
figure('Name', '40 - 90% 1RM Bland-Altman Plot (Left)')
tiledlayout(1,2);
nexttile;
scatter(left_average, left_difference, 'filled', 'MarkerFaceColor', [0 0 1]); hold on;
yline(left_upper, '--r');
textscatter(1.45, left_upper+0.01, "+1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
textscatter(1.45, left_upper-0.01, string(num2str(left_upper)), 'ColorData', [1 0 0])% Change first value accordingly
yline(left_average_difference_90, 'b');
textscatter(1.45, left_average_difference+0.01, "Mean", 'ColorData', [0 0 1]) % Change first value accordingly
textscatter(1.45, left_average_difference-0.01, string(num2str(left_average_difference)), 'ColorData', [0 0 1]) % Change first value accordingly
yline(left_lower, '--r');
textscatter(1.45, left_lower+0.01, "-1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
textscatter(1.45, left_lower-0.01, string(num2str(left_lower)), 'ColorData', [1 0 0]) % Change first value accordingly
title('40 - 90% 1RM Bland-Altman Plot (Left)'); xlabel('Mean of FLEX and VICON'); ylabel('Difference between FLEX and VICON'); hold off;

% Plot correlation graph
nexttile;
% figure('Name', '90% 1RM Correlation Plot (Left)')
scatter(all_left(:,1), all_left(:,2), 'filled');
lsline();
title('40 - 90% 1RM Correlation Plot (Left)'); xlabel('FLEX - Mean Concentric Velocity (m/s)'); ylabel('VICON - Mean Concentric Velocity (m/s)')

left_R = corrcoef(all_left);
sprintf('40 - 90%% Pearsons Correlation r (Left): %f', left_R(2))
%% Plot compile correlation table
Relative_Intensity = {'40%'; '60%';'70%';'80%';'90%'};
Left_Correlation = [left_R_40(2); left_R_60(2); left_R_70(2); left_R_80(2); left_R_90(2)];
Right_Correlation = [right_R_40(2); right_R_60(2); right_R_70(2); right_R_80(2); right_R_90(2)];
corr_table = table(Relative_Intensity, Left_Correlation, Right_Correlation);
%%
% alldata = array2table(alldata);
% alldata.Properties.VariableNames(1:end) = {'ParticipantID', '1RM', 'RepetitionNumber', 'FLEX_Mean_L', 'FLEX_Mean_R', ... 
%     'FLEX_Peak_L', 'FLEX_Peak_R', 'VICON_Mean_L', 'VICON_Mean_R', 'VICON_Peak_L','VICON_Peak_R'};
% lme = fitlme(alldata,'FLEX_Peak_L~VICON_Peak_L+(ParticipantID|1RM)');