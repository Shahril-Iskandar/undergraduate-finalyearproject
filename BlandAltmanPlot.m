clear
clc
close all
cd 'C:\Users\14000\Downloads\Year 4 FYP\Recruitment\Data Collection'
tic
%%

% Read dataset excel
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

% Combine mean data from both exercises
meanvelocitydata = vertcat(powerclean_meandata, snatch_meandata);
peakvelocitydata = vertcat(powerclean_peakdata, snatch_peakdata);
alldatacombined = vertcat(meanvelocitydata, peakvelocitydata);
%%
% Power Clean_Mean - Bland Altman
powerclean_meanvelocity_mean = zeros(6,3);
ptable = zeros(6,2);
figure('Name','Power Clean Mean - 1RM Bland-Altman Plot');
for i=1:length(intensity)
    data = powerclean_meandata(powerclean_meandata(:,2) == intensity(i),:); % Get values of specific intensity according to loop

    % Check normality
%     [H, pValue, W] = swtest(data(:,7),0.05);
%     if pValue <= 0.05
%         data = log(data);
%     end

    meandifferences = mean(data(:,7)); % Calculate the mean of the differences
    lowersd = meandifferences-1.96*std(data(:,7)); % Calculate lower LoA
    uppersd = meandifferences+1.96*std(data(:,7)); % Calculate upper LoA

%     mdl = fitlm(data(:,6),data(:,7)); % Linear regression
    [r,p] = corrcoef(data(:,6),data(:,7)); % Calculate r and p value of correlation
    sprintf('The r value for %d%% is %f', intensity(i), r(2)) 
    sprintf('The p value for %d%% is %f', intensity(i), p(2))
    ptable(i,1) = r(2); % Allocate r value
    ptable(i,2) = p(2); % Allocate p value
    if p(2) < 0.05
        sprintf('There is proportional bias for %d%% 1RM', intensity(i))
    else
        sprintf('There is no proportional bias for %d%% 1RM', intensity(i))
    end
  
    sprintf('The mean for %d%% is %f m/s', intensity(i), meandifferences)
    sprintf('Limits of agreement for %d%% is %f m/s to %f m/s', intensity(i), lowersd, uppersd)
    slope = data(:,4)\data(:,5);
    sprintf('Slope / Regression coefficient: %f', slope);
    powerclean_meanvelocity_mean(i,1) = meandifferences;
    powerclean_meanvelocity_mean(i,2) = lowersd;
    powerclean_meanvelocity_mean(i,3) = uppersd;
%     figure('Name',sprintf('%d%% 1RM Bland-Altman Plot', intensity(i)));
    subplot(2,3,i);
    scatter(data(:,6), data(:,7), 10, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410]);  % Plot mean velocity Bland Altman 
%     hold on;
    yline(uppersd, '--r');
    % textscatter(1.21, left_upper_40+0.01, "+1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
    % textscatter(1.21, left_upper_40-0.01, string(num2str(left_upper_40)), 'ColorData', [1 0 0])% Change first value accordingly
    yline(meandifferences, 'r');
    % textscatter(1.21, left_average_difference_40+0.01, "Mean", 'ColorData', [0 0 1]) % Change first value accordingly
    % textscatter(1.21, left_average_difference_40-0.01, string(num2str(left_average_difference_40)), 'ColorData', [0 0 1]) % Change first value accordingly
    yline(lowersd, '--r');
    % textscatter(1.21, left_lower_40+0.01, "-1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
    % textscatter(1.21, left_lower_40-0.01, string(num2str(left_lower_40)), 'ColorData', [1 0 0]) % Change first value accordingly
    % axis equal 
%     xlim([0.6 1.6]); ylim([-0.3 0.3]);
    L = lsline();
    L.Color = 'black';
    title(sprintf('%d%% 1RM', intensity(i))); xlabel('Mean of FLEX and VICON (m/s)'); ylabel('Difference (FLEX - VICON) (m/s)');
end
sgtitle('Power Clean - Mean Velocity (m/s)', 'FontSize', 13)
powerclean_meanvelocity_mean = [transpose(intensity) powerclean_meanvelocity_mean];
powerclean_meanvelocity_mean = array2table(powerclean_meanvelocity_mean, 'VariableNames', {'Intensity', 'Mean','LowerSD','UpperSD'});
disp('Power Clean Mean Velocity');
disp(powerclean_meanvelocity_mean);
ptable = [transpose(intensity) ptable];
ptable = array2table(ptable, 'VariableNames', {'Intensity', 'r-value', 'p-value'});
disp(ptable);

writetable(powerclean_meanvelocity_mean, 'AgreementLimit.xlsx','Sheet', 'Power Clean (Mean)')
%%
% Power Clean_Peak - Bland Altman
powerclean_peakvelocity_mean = zeros(6,3);
ptable = zeros(6,2);
figure('Name','Power Clean Peak - 1RM Bland-Altman Plot');
for i=1:length(intensity)
    data = powerclean_peakdata(powerclean_peakdata(:,2) == intensity(i),:); % Get values of specific intensity according to loop

    % Check normality
%     [H, pValue, W] = swtest(data(:,7),0.05);
%     if pValue <= 0.05
%         data = log(data);
%     end

    meandifferences = mean(data(:,7)); % Calculate the mean of the differences
    lowersd = meandifferences-1.96*std(data(:,7)); % Calculate lower LoA
    uppersd = meandifferences+1.96*std(data(:,7)); % Calculate upper LoA

%     mdl = fitlm(data(:,6),data(:,7)); % Linear regression
    [r,p] = corrcoef(data(:,6),data(:,7)); % Calculate r and p value of correlation
    sprintf('The r value for %d%% is %f', intensity(i), r(2)) 
    sprintf('The p value for %d%% is %f', intensity(i), p(2))
    ptable(i,1) = r(2); % Allocate r value
    ptable(i,2) = p(2); % Allocate p value
    if p(2) < 0.05
        sprintf('There is proportional bias for %d%% 1RM', intensity(i))
    else
        sprintf('There is no proportional bias for %d%% 1RM', intensity(i))
    end
  
    sprintf('The mean for %d%% is %f m/s', intensity(i), meandifferences)
    sprintf('Limits of agreement for %d%% is %f m/s to %f m/s', intensity(i), lowersd, uppersd)
    slope = data(:,4)\data(:,5);
    sprintf('Slope / Regression coefficient: %f', slope);
    powerclean_peakvelocity_mean(i,1) = meandifferences;
    powerclean_peakvelocity_mean(i,2) = lowersd;
    powerclean_peakvelocity_mean(i,3) = uppersd;
%     figure('Name',sprintf('%d%% 1RM Bland-Altman Plot', intensity(i)));
    subplot(2,3,i);
    scatter(data(:,6), data(:,7), 10, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410]);  % Plot mean velocity Bland Altman 
%     hold on;
    yline(uppersd, '--r');
    % textscatter(1.21, left_upper_40+0.01, "+1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
    % textscatter(1.21, left_upper_40-0.01, string(num2str(left_upper_40)), 'ColorData', [1 0 0])% Change first value accordingly
    yline(meandifferences, 'r');
    % textscatter(1.21, left_average_difference_40+0.01, "Mean", 'ColorData', [0 0 1]) % Change first value accordingly
    % textscatter(1.21, left_average_difference_40-0.01, string(num2str(left_average_difference_40)), 'ColorData', [0 0 1]) % Change first value accordingly
    yline(lowersd, '--r');
    % textscatter(1.21, left_lower_40+0.01, "-1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
    % textscatter(1.21, left_lower_40-0.01, string(num2str(left_lower_40)), 'ColorData', [1 0 0]) % Change first value accordingly
    % axis equal 
%     xlim([1.4 3.0]); ylim([-0.6 0.7])
    L = lsline();
    L.Color = 'black';
    title(sprintf('%d%% 1RM', intensity(i))); xlabel('Mean of FLEX and VICON (m/s)'); ylabel('Difference (FLEX - VICON) (m/s)');
end
sgtitle('Power Clean - Peak Velocity (m/s)', 'FontSize', 13)
powerclean_peakvelocity_mean = [transpose(intensity) powerclean_peakvelocity_mean];
powerclean_peakvelocity_mean = array2table(powerclean_peakvelocity_mean, 'VariableNames', {'Intensity', 'Mean','LowerSD','UpperSD'});
disp('Power Clean Peak Velocity');
disp(powerclean_peakvelocity_mean);
ptable = [transpose(intensity) ptable];
ptable = array2table(ptable, 'VariableNames', {'Intensity', 'r-value', 'p-value'});
disp(ptable);

writetable(powerclean_peakvelocity_mean, 'AgreementLimit.xlsx','Sheet', 'Power Clean (Peak)')
%%
% Snatch_Mean- Bland Altman
snatch_meanvelocity_mean = zeros(6,3);
ptable = zeros(6,2);
figure('Name','Snatch Mean - 1RM Bland-Altman Plot');
for i=1:length(intensity)
    data = snatch_meandata(snatch_meandata(:,2) == intensity(i),:); % Get values of specific intensity according to loop

    % Check normality
%     [H, pValue, W] = swtest(data(:,7),0.05);
%     if pValue <= 0.05
%         data = log(data);
%     end

    meandifferences = mean(data(:,7)); % Calculate the mean of the differences
    lowersd = meandifferences-1.96*std(data(:,7)); % Calculate lower LoA
    uppersd = meandifferences+1.96*std(data(:,7)); % Calculate upper LoA

%     mdl = fitlm(data(:,6),data(:,7)); % Linear regression
    [r,p] = corrcoef(data(:,6),data(:,7)); % Calculate r and p value of correlation
    sprintf('The r value for %d%% is %f', intensity(i), r(2)) 
    sprintf('The p value for %d%% is %f', intensity(i), p(2))
    ptable(i,1) = r(2); % Allocate r value
    ptable(i,2) = p(2); % Allocate p value
    if p(2) < 0.05
        sprintf('There is proportional bias for %d%% 1RM', intensity(i))
    else
        sprintf('There is no proportional bias for %d%% 1RM', intensity(i))
    end
  
    sprintf('The mean for %d%% is %f m/s', intensity(i), meandifferences)
    sprintf('Limits of agreement for %d%% is %f m/s to %f m/s', intensity(i), lowersd, uppersd)
    slope = data(:,4)\data(:,5);
    sprintf('Slope / Regression coefficient: %f', slope);
    snatch_meanvelocity_mean(i,1) = meandifferences;
    snatch_meanvelocity_mean(i,2) = lowersd;
    snatch_meanvelocity_mean(i,3) = uppersd;
%     figure('Name',sprintf('%d%% 1RM Bland-Altman Plot', intensity(i)));
    subplot(2,3,i);
    scatter(data(:,6), data(:,7),  10, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410]);  % Plot mean velocity Bland Altman 
%     hold on;
    yline(uppersd, '--r');
    % textscatter(1.21, left_upper_40+0.01, "+1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
    % textscatter(1.21, left_upper_40-0.01, string(num2str(left_upper_40)), 'ColorData', [1 0 0])% Change first value accordingly
    yline(meandifferences, 'r');
    % textscatter(1.21, left_average_difference_40+0.01, "Mean", 'ColorData', [0 0 1]) % Change first value accordingly
    % textscatter(1.21, left_average_difference_40-0.01, string(num2str(left_average_difference_40)), 'ColorData', [0 0 1]) % Change first value accordingly
    yline(lowersd, '--r');
    % textscatter(1.21, left_lower_40+0.01, "-1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
    % textscatter(1.21, left_lower_40-0.01, string(num2str(left_lower_40)), 'ColorData', [1 0 0]) % Change first value accordingly
    % axis equal 
%     xlim([0.7 1.6]); ylim([-0.4 0.3])
    L = lsline();
    L.Color = 'black';
    title(sprintf('%d%% 1RM', intensity(i))); xlabel('Mean of FLEX and VICON (m/s)'); ylabel('Difference (FLEX - VICON) (m/s)');
end
sgtitle('Snatch - Mean Velocity (m/s)', 'FontSize', 13)
snatch_meanvelocity_mean = [transpose(intensity) snatch_meanvelocity_mean];
snatch_meanvelocity_mean = array2table(snatch_meanvelocity_mean, 'VariableNames', {'Intensity', 'Mean','LowerSD','UpperSD'});
disp('Snatch Mean Velocity');
disp(snatch_meanvelocity_mean);
ptable = [transpose(intensity) ptable];
ptable = array2table(ptable, 'VariableNames', {'Intensity', 'r-value', 'p-value'});
disp(ptable);

writetable(snatch_meanvelocity_mean, 'AgreementLimit.xlsx','Sheet', 'Snatch (Mean)')
%%
% Snatch_Peak - Bland Altman
snatch_peakvelocity_mean = zeros(6,3);
ptable = zeros(6,2);
figure('Name','Snatch Peak - 1RM Bland-Altman Plot');
for i=1:length(intensity)
    data = snatch_peakdata(snatch_peakdata(:,2) == intensity(i),:); % Get values of specific intensity according to loop

    % Check normality
%     [H, pValue, W] = swtest(data(:,7),0.05);
%     if pValue <= 0.05
%         data = log(data);
%     end

    meandifferences = mean(data(:,7)); % Calculate the mean of the differences
    lowersd = meandifferences-1.96*std(data(:,7)); % Calculate lower LoA
    uppersd = meandifferences+1.96*std(data(:,7)); % Calculate upper LoA

%     mdl = fitlm(data(:,6),data(:,7)); % Linear regression
    [r,p] = corrcoef(data(:,6),data(:,7)); % Calculate r and p value of correlation
    sprintf('The r value for %d%% is %f', intensity(i), r(2)) 
    sprintf('The p value for %d%% is %f', intensity(i), p(2))
    ptable(i,1) = r(2); % Allocate r value
    ptable(i,2) = p(2); % Allocate p value
    if p(2) < 0.05
        sprintf('There is proportional bias for %d%% 1RM', intensity(i))
    else
        sprintf('There is no proportional bias for %d%% 1RM', intensity(i))
    end
  
    sprintf('The mean for %d%% is %f m/s', intensity(i), meandifferences)
    sprintf('Limits of agreement for %d%% is %f m/s to %f m/s', intensity(i), lowersd, uppersd)
    slope = data(:,4)\data(:,5);
    sprintf('Slope / Regression coefficient: %f', slope);
    snatch_peakvelocity_mean(i,1) = meandifferences;
    snatch_peakvelocity_mean(i,2) = lowersd;
    snatch_peakvelocity_mean(i,3) = uppersd;
%     figure('Name',sprintf('%d%% 1RM Bland-Altman Plot', intensity(i)));
    subplot(2,3,i);
    scatter(data(:,6), data(:,7), 10, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410]);  % Plot mean velocity Bland Altman 
%     hold on;
    yline(uppersd, '--r');
    % textscatter(1.21, left_upper_40+0.01, "+1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
    % textscatter(1.21, left_upper_40-0.01, string(num2str(left_upper_40)), 'ColorData', [1 0 0])% Change first value accordingly
    yline(meandifferences, 'r');
    % textscatter(1.21, left_average_difference_40+0.01, "Mean", 'ColorData', [0 0 1]) % Change first value accordingly
    % textscatter(1.21, left_average_difference_40-0.01, string(num2str(left_average_difference_40)), 'ColorData', [0 0 1]) % Change first value accordingly
    yline(lowersd, '--r');
    % textscatter(1.21, left_lower_40+0.01, "-1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
    % textscatter(1.21, left_lower_40-0.01, string(num2str(left_lower_40)), 'ColorData', [1 0 0]) % Change first value accordingly
    % axis equal 
%     xlim([1.5 3.2]); ylim([-1.5 0.5])
    L = lsline();
    L.Color = 'black';
    title(sprintf('%d%% 1RM', intensity(i))); xlabel('Mean of FLEX and VICON (m/s)'); ylabel('Difference (FLEX - VICON) (m/s)');
end
sgtitle('Snatch - Peak Velocity (m/s)', 'FontSize', 13)
snatch_peakvelocity_mean = [transpose(intensity) snatch_peakvelocity_mean];
snatch_peakvelocity_mean = array2table(snatch_peakvelocity_mean, 'VariableNames', {'Intensity', 'Mean','LowerSD','UpperSD'});
disp('Snatch - Peak Velocity');
disp(snatch_peakvelocity_mean);
ptable = [transpose(intensity) ptable];
ptable = array2table(ptable, 'VariableNames', {'Intensity', 'r-value', 'p-value'});
disp(ptable);

writetable(snatch_peakvelocity_mean, 'AgreementLimit.xlsx','Sheet', 'Snatch (Peak)')
%% Plot all intensity

% Power Clean mean all intensity
allintensity = zeros(4,3);
ptable = zeros(1,2);
figure('Name','Power Clean Mean (All %1RM) Bland-Altman Plot');
data = powerclean_meandata;
meandifferences = mean(data(:,7)); % Calculate the mean of the differences
lowersd = meandifferences-1.96*std(data(:,7)); % Calculate lower LoA
uppersd = meandifferences+1.96*std(data(:,7)); % Calculate upper LoA
%     mdl = fitlm(data(:,6),data(:,7)); % Linear regression
[r,p] = corrcoef(data(:,6),data(:,7)); % Calculate r and p value of correlation
sprintf('The r value for %d%% is %f', intensity(i), r(2)) 
sprintf('The p value for %d%% is %f', intensity(i), p(2))
ptable(1,1) = r(2); % Allocate r value
ptable(1,2) = p(2); % Allocate p value
if p(2) < 0.05
    disp('There is proportional bias for Power Clean - Mean Velocity (Overall)');
else
    disp('There is no proportional bias for Power Clean - Mean Velocity (Overall)');
end
sprintf('The mean for Power Clean Mean (All %%1RM): %f m/s', meandifferences);
sprintf('Limits of agreement for Power Clean Mean (All %%1RM): %f m/s to %f m/s', lowersd, uppersd);
allintensity(1,1) = meandifferences;
allintensity(1,2) = lowersd;
allintensity(1,3) = uppersd;
slope = data(:,4)\data(:,5);
sprintf('Slope / Regression coefficient: %f', slope);
scatter(data(:,6), data(:,7), 15, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410]);  % Plot mean velocity Bland Altman
% gscatter(data(:,6), data(:,7), data(:,2), 'krbmgc')
yline(uppersd, '--r', 'LineWidth', 2, 'HandleVisibility','off');
% textscatter(1.21, left_upper_40+0.01, "+1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
% textscatter(1.21, left_upper_40-0.01, string(num2str(left_upper_40)), 'ColorData', [1 0 0])% Change first value accordingly
yline(meandifferences, 'r', 'LineWidth', 2, 'HandleVisibility','off');
% textscatter(1.21, left_average_difference_40+0.01, "Mean", 'ColorData', [0 0 1]) % Change first value accordingly
% textscatter(1.21, left_average_difference_40-0.01, string(num2str(left_average_difference_40)), 'ColorData', [0 0 1]) % Change first value accordingly
yline(lowersd, '--r', 'LineWidth', 2, 'HandleVisibility','off');
% textscatter(1.21, left_lower_40+0.01, "-1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
% textscatter(1.21, left_lower_40-0.01, string(num2str(left_lower_40)), 'ColorData', [1 0 0]) % Change first value accordingly
L = lsline();
L.Color = 'black';
title('Power Clean - Mean Velocity (All %1RM)'); xlabel('Mean of FLEX and VICON (m/s)'); ylabel('Difference (FLEX - VICON) (m/s)');
ax = gca;
ax.FontSize = 16;
ptable = array2table(ptable, 'VariableNames', {'r-value', 'p-value'});
disp(ptable);

% tscore = tinv(.95,size(data,2));
% standarderror = std(data(:,7))/sqrt(size(data,2));
% upperconfidence = meandifferences + tscore * standarderror;
% lowerconfidence = meandifferences - tscore * standarderror;
sprintf('Mean difference is: %f', meandifferences)
sprintf('95CI is: %f to %f', lowersd, uppersd)

% Power Clean peak all intensity
figure('Name','Power Clean Peak (All %1RM) Bland-Altman Plot');
ptable = zeros(1,2);
data = powerclean_peakdata;
meandifferences = mean(data(:,7)); % Calculate the mean of the differences
lowersd = meandifferences-1.96*std(data(:,7)); % Calculate lower LoA
uppersd = meandifferences+1.96*std(data(:,7)); % Calculate upper LoA
%     mdl = fitlm(data(:,6),data(:,7)); % Linear regression
[r,p] = corrcoef(data(:,6),data(:,7)); % Calculate r and p value of correlation
sprintf('The r value for %d%% is %f', intensity(i), r(2)) 
sprintf('The p value for %d%% is %f', intensity(i), p(2))
ptable(1,1) = r(2); % Allocate r value
ptable(1,2) = p(2); % Allocate p value
if p(2) < 0.05
    disp('There is proportional bias for Power Clean - Peak Velocity (Overall)');
else
    disp('There is no proportional bias for Power Clean - Peak Velocity (Overall)');
end
sprintf('The mean for Power Clean Mean (All %%1RM): %f m/s', meandifferences);
sprintf('Limits of agreement for Power Clean Mean (All %%1RM): %f m/s to %f m/s', lowersd, uppersd);
allintensity(2,1) = meandifferences;
allintensity(2,2) = lowersd;
allintensity(2,3) = uppersd;
slope = data(:,4)\data(:,5);
sprintf('Slope / Regression coefficient: %f', slope);
scatter(data(:,6), data(:,7), 15, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410]);  % Plot mean velocity Bland Altman 
yline(uppersd, '--r','LineWidth', 2, 'HandleVisibility','off');
% textscatter(1.21, left_upper_40+0.01, "+1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
% textscatter(1.21, left_upper_40-0.01, string(num2str(left_upper_40)), 'ColorData', [1 0 0])% Change first value accordingly
yline(meandifferences, 'r','LineWidth', 2, 'HandleVisibility','off');
% textscatter(1.21, left_average_difference_40+0.01, "Mean", 'ColorData', [0 0 1]) % Change first value accordingly
% textscatter(1.21, left_average_difference_40-0.01, string(num2str(left_average_difference_40)), 'ColorData', [0 0 1]) % Change first value accordingly
yline(lowersd, '--r','LineWidth', 2, 'HandleVisibility','off');
% textscatter(1.21, left_lower_40+0.01, "-1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
% textscatter(1.21, left_lower_40-0.01, string(num2str(left_lower_40)), 'ColorData', [1 0 0]) % Change first value accordingly
L = lsline();
L.Color = 'black';
title('Power Clean - Peak Velocity (All %1RM)'); xlabel('Mean of FLEX and VICON (m/s)'); ylabel('Difference (FLEX - VICON) (m/s)');
ax = gca;
ax.FontSize = 16;
ptable = array2table(ptable, 'VariableNames', {'r-value', 'p-value'});
disp(ptable);
sprintf('Mean difference is: %f', meandifferences)
sprintf('95CI is: %f to %f', lowersd, uppersd)

% Snatch mean all intensity
figure('Name','Snatch Mean (All %1RM) Bland-Altman Plot');
ptable = zeros(1,2);
data = snatch_meandata;
meandifferences = mean(data(:,7)); % Calculate the mean of the differences
lowersd = meandifferences-1.96*std(data(:,7)); % Calculate lower LoA
uppersd = meandifferences+1.96*std(data(:,7)); % Calculate upper LoA
%     mdl = fitlm(data(:,6),data(:,7)); % Linear regression
[r,p] = corrcoef(data(:,6),data(:,7)); % Calculate r and p value of correlation
sprintf('The r value for %d%% is %f', intensity(i), r(2)) 
sprintf('The p value for %d%% is %f', intensity(i), p(2))
ptable(1,1) = r(2); % Allocate r value
ptable(1,2) = p(2); % Allocate p value
if p(2) < 0.05
    disp('There is proportional bias for Snatch - Mean Velocity (Overall)');
else
    disp('There is no proportional bias for Snatch - Mean Velocity (Overall)');
end
sprintf('The mean for Snatch Mean (All %%1RM): %f m/s', meandifferences);
sprintf('Limits of agreement for Snatch Mean (All %%1RM): %f m/s to %f m/s', lowersd, uppersd);
allintensity(3,1) = meandifferences;
allintensity(3,2) = lowersd;
allintensity(3,3) = uppersd;
slope = data(:,4)\data(:,5);
sprintf('Slope / Regression coefficient: %f', slope);
scatter(data(:,6), data(:,7),  15, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410]);  % Plot mean velocity Bland Altman 
yline(uppersd, '--r','LineWidth', 2, 'HandleVisibility','off');
% textscatter(1.21, left_upper_40+0.01, "+1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
% textscatter(1.21, left_upper_40-0.01, string(num2str(left_upper_40)), 'ColorData', [1 0 0])% Change first value accordingly
yline(meandifferences, 'r','LineWidth', 2, 'HandleVisibility','off');
% textscatter(1.21, left_average_difference_40+0.01, "Mean", 'ColorData', [0 0 1]) % Change first value accordingly
% textscatter(1.21, left_average_difference_40-0.01, string(num2str(left_average_difference_40)), 'ColorData', [0 0 1]) % Change first value accordingly
yline(lowersd, '--r','LineWidth', 2, 'HandleVisibility','off');
% textscatter(1.21, left_lower_40+0.01, "-1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
% textscatter(1.21, left_lower_40-0.01, string(num2str(left_lower_40)), 'ColorData', [1 0 0]) % Change first value accordingly
L = lsline();
L.Color = 'black';
title('Snatch - Mean Velocity (All %1RM)'); xlabel('Mean of FLEX and VICON (m/s)'); ylabel('Difference (FLEX - VICON) (m/s)');
ax = gca;
ax.FontSize = 16;
ptable = array2table(ptable, 'VariableNames', {'r-value', 'p-value'});
disp(ptable);
sprintf('Mean difference is: %f', meandifferences)
sprintf('95CI is: %f to %f', lowersd, uppersd)

% Snatch peak all intensity
figure('Name','Snatch Peak (All %1RM) Bland-Altman Plot');
ptable = zeros(1,2);
data = snatch_peakdata;
meandifferences = mean(data(:,7)); % Calculate the mean of the differences
lowersd = meandifferences-1.96*std(data(:,7)); % Calculate lower LoA
uppersd = meandifferences+1.96*std(data(:,7)); % Calculate upper LoA
%     mdl = fitlm(data(:,6),data(:,7)); % Linear regression
[r,p] = corrcoef(data(:,6),data(:,7)); % Calculate r and p value of correlation
sprintf('The r value for %d%% is %f', intensity(i), r(2)) 
sprintf('The p value for %d%% is %f', intensity(i), p(2))
ptable(1,1) = r(2); % Allocate r value
ptable(1,2) = p(2); % Allocate p value
if p(2) < 0.05
    disp('There is proportional bias for Snatch - Peak Velocity (Overall)');
else
    disp('There is no proportional bias for Snatch - Peak DVelocity (Overall)');
end
sprintf('The mean for Snatch Peak (All %%1RM): %f m/s', meandifferences);
sprintf('Limits of agreement for Snatch Peak (All %%1RM): %f m/s to %f m/s', lowersd, uppersd);
allintensity(4,1) = meandifferences;
allintensity(4,2) = lowersd;
allintensity(4,3) = uppersd;
slope = data(:,4)\data(:,5);
sprintf('Slope / Regression coefficient: %f', slope);
scatter(data(:,6), data(:,7),  15, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410]);  % Plot mean velocity Bland Altman 
yline(uppersd, '--r','LineWidth', 2, 'HandleVisibility','off');
% textscatter(1.21, left_upper_40+0.01, "+1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
% textscatter(1.21, left_upper_40-0.01, string(num2str(left_upper_40)), 'ColorData', [1 0 0])% Change first value accordingly
yline(meandifferences, 'r','LineWidth', 2, 'HandleVisibility','off');
% textscatter(1.21, left_average_difference_40+0.01, "Mean", 'ColorData', [0 0 1]) % Change first value accordingly
% textscatter(1.21, left_average_difference_40-0.01, string(num2str(left_average_difference_40)), 'ColorData', [0 0 1]) % Change first value accordingly
yline(lowersd, '--r','LineWidth', 2, 'HandleVisibility','off');
% textscatter(1.21, left_lower_40+0.01, "-1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
% textscatter(1.21, left_lower_40-0.01, string(num2str(left_lower_40)), 'ColorData', [1 0 0]) % Change first value accordingly
L = lsline();
L.Color = 'black';
title('Snatch - Peak Velocity (All %1RM)'); xlabel('Mean of FLEX and VICON (m/s)'); ylabel('Difference (FLEX - VICON) (m/s)');
ax = gca;
ax.FontSize = 16;
sprintf('Mean difference is: %f', meandifferences)
sprintf('95CI is: %f to %f', lowersd, uppersd)

condition = ["Power Clean (Mean)"; "Power Clean (Peak)"; "Snatch (Mean)"; "Snatch (Peak)"];
allintensity = [condition allintensity];
allintensity = array2table(allintensity, 'VariableNames', {'Intensity', 'Mean','LowerSD','UpperSD'});
disp('All Intensities Combined in each exercise');
disp(allintensity);
ptable = array2table(ptable, 'VariableNames', {'r-value', 'p-value'});
disp(ptable);
%% Power clean and snatch mea data combined, different colors

% Mean concentric velocity
figure('Name','Mean velocity all exercises (All %1RM) Bland-Altman Plot');
ptable = zeros(1,2);
data = meanvelocitydata;
meandifferences = mean(data(:,7)); % Calculate the mean of the differences
lowersd = meandifferences-1.96*std(data(:,7)); % Calculate lower LoA
uppersd = meandifferences+1.96*std(data(:,7)); % Calculate upper LoA
h1 = scatter(data(:,6), data(:,7),  13, 'filled');  % Plot mean velocity Bland Altman
hold on;
L = lsline();
L.Color = 'black';
L.LineWidth = 3;
data = powerclean_meandata;
h2 = scatter(data(:,6), data(:,7),  20, 'filled', 'MarkerFaceColor', "#FF0000");  % Plot mean velocity Bland Altman
data = snatch_meandata;
h3 = scatter(data(:,6), data(:,7),  20, 'filled', 'MarkerFaceColor', "#0000FF");  % Plot mean velocity Bland Altman
y1 = yline(uppersd, '--k', num2str(round(uppersd,2)),'HandleVisibility','off', 'LineWidth', 3);
y2 = yline(meandifferences, 'k', num2str(round(meandifferences,2)),'HandleVisibility','off','LineWidth', 3);
y3 = yline(lowersd, '--k', num2str(round(lowersd,2)),'HandleVisibility','off','LineWidth', 3);
title('Mean Velocity from both exercises (All %1RM)', 'FontSize', 13); xlabel('Mean of FLEX and VICON (m/s)'); ylabel('Difference (FLEX - VICON) (m/s)');
y1.FontSize = 16;
y2.FontSize = 16;
y3.FontSize = 16;
legend([h2 h3], {'Power Clean', 'Snatch'});
set(gca,'FontSize',22);

% Peak concentric velocity 
figure('Name','Mean velocity all exercises (All %1RM) Bland-Altman Plot');
ptable = zeros(1,2);
data = peakvelocitydata;
meandifferences = mean(data(:,7)); % Calculate the mean of the differences
lowersd = meandifferences-1.96*std(data(:,7)); % Calculate lower LoA
uppersd = meandifferences+1.96*std(data(:,7)); % Calculate upper LoA
h1 = scatter(data(:,6), data(:,7),  20, 'filled');  % Plot mean velocity Bland Altman
hold on;
L = lsline();
L.Color = 'green';
L.LineWidth = 3;
data = powerclean_peakdata;
h2 = scatter(data(:,6), data(:,7),  20, 'filled', 'MarkerFaceColor', "#FF0000");  % Plot mean velocity Bland Altman
data = snatch_peakdata;
h3 = scatter(data(:,6), data(:,7),  20, 'filled', 'MarkerFaceColor', "#0000FF");  % Plot mean velocity Bland Altman
y1 = yline(uppersd, '--k', num2str(round(uppersd,2)),'HandleVisibility','off', 'LineWidth', 3);
y2 = yline(meandifferences, 'k', num2str(round(meandifferences,2)),'HandleVisibility','off','LineWidth', 3);
y3 = yline(lowersd, '--k', num2str(round(lowersd,2)),'HandleVisibility','off','LineWidth', 3);
y1.FontSize = 16;
y2.FontSize = 16;
y3.FontSize = 16;
title('Peak Velocity (All %1RM)', 'FontSize', 13); xlabel('Mean of FLEX and VICON (m/s)'); ylabel('Difference (FLEX - VICON) (m/s)');
legend([h2 h3], {'Power Clean', 'Snatch'});
set(gca,'FontSize',22);
%% Mean velocity data combined from power clean and snatch

% subplot(3,1,1)
% Mean velocity
figure('Name','Mean velocity all exercises (All %1RM) Bland-Altman Plot');
ptable = zeros(1,2);
data = meanvelocitydata;
meandifferences = mean(data(:,7)); % Calculate the mean of the differences
lowersd = meandifferences-1.96*std(data(:,7)); % Calculate lower LoA
uppersd = meandifferences+1.96*std(data(:,7)); % Calculate upper LoA
%     mdl = fitlm(data(:,6),data(:,7)); % Linear regression
[r,p] = corrcoef(data(:,6),data(:,7)); % Calculate r and p value of correlation
sprintf('The r value for %d%% is %f', intensity(i), r(2)) 
sprintf('The p value for %d%% is %f', intensity(i), p(2)) 
ptable(1,1) = r(2); % Allocate r value
ptable(1,2) = p(2); % Allocate p value
if p(2) < 0.05
    disp('There is proportional bias for all exercises (Overall)');
else
    disp('There is no proportional bias for all exercises  Overall)');
end
sprintf('The mean for all exercises (All %%1RM): %f m/s', meandifferences);
sprintf('Limits of agreement for all exercises (All %%1RM): %f m/s to %f m/s', lowersd, uppersd);
% allintensity(4,1) = meandifferences;
% allintensity(4,2) = lowersd;
% allintensity(4,3) = uppersd;
slope = data(:,4)\data(:,5);
sprintf('Slope / Regression coefficient: %f', slope);
% gscatter(data(:,6), data(:,7), data(:,2), 'krbmgc');
scatter(data(:,6), data(:,7),  10, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410]);  % Plot mean velocity Bland Altman 
% yline([uppersd meandifferences lowersd], '--r', {'uppersd', 'mean', 'lowersd'},'HandleVisibility','off')
yline(uppersd, '--r','HandleVisibility','off');
% textscatter(1.21, left_upper_40+0.01, "+1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
% textscatter(1.21, left_upper_40-0.01, string(num2str(left_upper_40)), 'ColorData', [1 0 0])% Change first value accordingly
yline(meandifferences, 'r','HandleVisibility','off');
% textscatter(1.21, left_average_difference_40+0.01, "Mean", 'ColorData', [0 0 1]) % Change first value accordingly
% textscatter(1.21, left_average_difference_40-0.01, string(num2str(left_average_difference_40)), 'ColorData', [0 0 1]) % Change first value accordingly
yline(lowersd, '--r' ,'HandleVisibility','off');
% textscatter(1.21, left_lower_40+0.01, "-1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
% textscatter(1.21, left_lower_40-0.01, string(num2str(left_lower_40)), 'ColorData', [1 0 0]) % Change first value accordingly
L = lsline();
L.Color = 'black';
title('Mean Velocity from both exercises (All %1RM)'); xlabel('Mean of FLEX and VICON (m/s)'); ylabel('Difference (FLEX - VICON) (m/s)');

% condition = ["Power Clean (Mean)"; "Power Clean (Peak)"; "Snatch (Mean)"; "Snatch (Peak)"];
% allintensity = [condition allintensity];
% allintensity = array2table(allintensity, 'VariableNames', {'Intensity', 'Mean','LowerSD','UpperSD'});
% disp('All Intensities Combined in each exercise');
% disp(allintensity);
% ptable = array2table(ptable, 'VariableNames', {'r-value', 'p-value'});
% disp(ptable);

% subplot(3,1,2)
% Peak velocity
figure('Name','Peak velocity all exercises (All %1RM) Bland-Altman Plot');
ptable = zeros(1,2);
data = peakvelocitydata;
meandifferences = mean(data(:,7)); % Calculate the mean of the differences
lowersd = meandifferences-1.96*std(data(:,7)); % Calculate lower LoA
uppersd = meandifferences+1.96*std(data(:,7)); % Calculate upper LoA
%     mdl = fitlm(data(:,6),data(:,7)); % Linear regression
[r,p] = corrcoef(data(:,6),data(:,7)); % Calculate r and p value of correlation
sprintf('The r value for %d%% is %f', intensity(i), r(2)) 
sprintf('The p value for %d%% is %f', intensity(i), p(2))
ptable(1,1) = r(2); % Allocate r value
ptable(1,2) = p(2); % Allocate p value
if p(2) < 0.05
    disp('There is proportional bias for all exercises - Peak (Overall)');
else
    disp('There is no proportional bias for all exercises - Peak Overall)');
end
sprintf('The mean for all exercises Peak (All %%1RM): %f m/s', meandifferences);
sprintf('Limits of agreement for all exercises Peak (All %%1RM): %f m/s to %f m/s', lowersd, uppersd);
% allintensity(4,1) = meandifferences;
% allintensity(4,2) = lowersd;
% allintensity(4,3) = uppersd;
slope = data(:,4)\data(:,5);
sprintf('Slope / Regression coefficient: %f', slope);
scatter(data(:,6), data(:,7),  10, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410]);  % Plot mean velocity Bland Altman 
yline(uppersd, '--r');
% textscatter(1.21, left_upper_40+0.01, "+1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
% textscatter(1.21, left_upper_40-0.01, string(num2str(left_upper_40)), 'ColorData', [1 0 0])% Change first value accordingly
yline(meandifferences, 'b');
% textscatter(1.21, left_average_difference_40+0.01, "Mean", 'ColorData', [0 0 1]) % Change first value accordingly
% textscatter(1.21, left_average_difference_40-0.01, string(num2str(left_average_difference_40)), 'ColorData', [0 0 1]) % Change first value accordingly
yline(lowersd, '--r');
% textscatter(1.21, left_lower_40+0.01, "-1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
% textscatter(1.21, left_lower_40-0.01, string(num2str(left_lower_40)), 'ColorData', [1 0 0]) % Change first value accordingly
L = lsline();
L.Color = 'black';
title('Peak Velocity from both exercises (All %1RM)'); xlabel('Mean of FLEX and VICON (m/s)'); ylabel('Difference (FLEX - VICON) (m/s)');

% condition = ["Power Clean (Mean)"; "Power Clean (Peak)"; "Snatch (Mean)"; "Snatch (Peak)"];
% allintensity = [condition allintensity];
% allintensity = array2table(allintensity, 'VariableNames', {'Intensity', 'Mean','LowerSD','UpperSD'});
% disp('All Intensities Combined in each exercise');
% disp(allintensity);
% ptable = array2table(ptable, 'VariableNames', {'r-value', 'p-value'});
% disp(ptable);

% All data mean and peak combined
% subplot(3,1,3)
figure('Name','All velocity all exercises (All %1RM) Bland-Altman Plot');
ptable = zeros(1,2);
data = alldatacombined;
meandifferences = mean(data(:,7)); % Calculate the mean of the differences
lowersd = meandifferences-1.96*std(data(:,7)); % Calculate lower LoA
uppersd = meandifferences+1.96*std(data(:,7)); % Calculate upper LoA
%     mdl = fitlm(data(:,6),data(:,7)); % Linear regression
[r,p] = corrcoef(data(:,6),data(:,7)); % Calculate r and p value of correlation
sprintf('The r value for %d%% is %f', intensity(i), r(2)) 
sprintf('The p value for %d%% is %f', intensity(i), p(2))
ptable(1,1) = r(2); % Allocate r value
ptable(1,2) = p(2); % Allocate p value
if p(2) < 0.05
    disp('There is proportional bias for all exercises - all (Overall)');
else
    disp('There is no proportional bias for all exercises - all Overall)');
end
sprintf('The mean for all exercises all (All %%1RM): %f m/s', meandifferences);
sprintf('Limits of agreement for all exercises all (All %%1RM): %f m/s to %f m/s', lowersd, uppersd);
% allintensity(4,1) = meandifferences;
% allintensity(4,2) = lowersd;
% allintensity(4,3) = uppersd;
slope = data(:,4)\data(:,5);
sprintf('Slope / Regression coefficient: %f', slope);
scatter(data(:,6), data(:,7),  10, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410]);  % Plot mean velocity Bland Altman 
yline(uppersd, '--r');
% textscatter(1.21, left_upper_40+0.01, "+1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
% textscatter(1.21, left_upper_40-0.01, string(num2str(left_upper_40)), 'ColorData', [1 0 0])% Change first value accordingly
yline(meandifferences, 'b');
% textscatter(1.21, left_average_difference_40+0.01, "Mean", 'ColorData', [0 0 1]) % Change first value accordingly
% textscatter(1.21, left_average_difference_40-0.01, string(num2str(left_average_difference_40)), 'ColorData', [0 0 1]) % Change first value accordingly
yline(lowersd, '--r');
% textscatter(1.21, left_lower_40+0.01, "-1.96 SD", 'ColorData', [1 0 0]) % Change first value accordingly
% textscatter(1.21, left_lower_40-0.01, string(num2str(left_lower_40)), 'ColorData', [1 0 0]) % Change first value accordingly
L = lsline();
L.Color = 'black';
title('All velocity from both exercises (All %1RM)'); xlabel('Mean of FLEX and VICON (m/s)'); ylabel('Difference (FLEX - VICON) (m/s)');

% condition = ["Power Clean (Mean)"; "Power Clean (Peak)"; "Snatch (Mean)"; "Snatch (Peak)"];
% allintensity = [condition allintensity];
% allintensity = array2table(allintensity, 'VariableNames', {'Intensity', 'Mean','LowerSD','UpperSD'});
% disp('All Intensities Combined in each exercise');
% disp(allintensity);
% ptable = array2table(ptable, 'VariableNames', {'r-value', 'p-value'});
% disp(ptable);


%% Correlation plot for all intensity

[powercleancorrelation_mean powercleancorrelation_mean_r]= corrcoef(powerclean_meandata(:,4), powerclean_meandata(:,5));
[powercleancorrelation_peak powercleancorrelation_peak_r] = corrcoef(powerclean_peakdata(:,4), powerclean_peakdata(:,5));
[snatchcorrelation_mean snatchcorrelation_mean_r] = corrcoef(snatch_meandata(:,4), snatch_meandata(:,5));
[snatchcorrelation_peak snatchcorrelation_peak_r] = corrcoef(snatch_peakdata(:,4), snatch_peakdata(:,5));

subplot(2,2,1)
% Power Clean Mean Data
scatter(powerclean_meandata(:,4), powerclean_meandata(:,5), 10, 'filled');
L = lsline();
L.Color = 'black';
title('Power Clean - Mean Velocity (All % 1RM)'); xlabel('FLEX - Mean Concentric Velocity (m/s)'); ylabel('VICON - Mean Concentric Velocity (m/s)');
text(1.4, 0.9, ['r = ' num2str(round(powercleancorrelation_mean(2),2))])

subplot(2,2,2)
% Power Clean Peak Data
scatter(powerclean_peakdata(:,4), powerclean_peakdata(:,5), 10, 'filled');
L = lsline();
L.Color = 'black';
title('Power Clean - Peak Velocity (All % 1RM)'); xlabel('FLEX - Mean Concentric Velocity (m/s)'); ylabel('VICON - Mean Concentric Velocity (m/s)');
text(3, 1.8, ['r = ' num2str(round(powercleancorrelation_peak(2),2))])

subplot(2,2,3)
% Snatch Mean Data
scatter(snatch_meandata(:,4), snatch_meandata(:,5), 10, 'filled');
L = lsline();
L.Color = 'black';
title('Snatch - Mean Velocity (All % 1RM)'); xlabel('FLEX - Mean Concentric Velocity (m/s)'); ylabel('VICON - Mean Concentric Velocity (m/s)');
text(1.4, 0.9, ['r = ' num2str(round(snatchcorrelation_mean(2),2))])

subplot(2,2,4)
% Power Clean Peak Data
scatter(snatch_peakdata(:,4), snatch_peakdata(:,5), 10, 'filled');
L = lsline();
L.Color = 'black';
title('Snatch - Peak Velocity (All % 1RM)'); xlabel('FLEX - Mean Concentric Velocity (m/s)'); ylabel('VICON - Mean Concentric Velocity (m/s)');
text(3, 1.8, ['r = ' num2str(round(snatchcorrelation_peak(2),2))])

correlation = [powercleancorrelation_mean(2) powercleancorrelation_peak(2) snatchcorrelation_mean(2) snatchcorrelation_peak(2)];
correlation = array2table(correlation, 'VariableNames', {'Power Clean - Mean Velocity', 'Power Clean - Peak Velocity', 'Snatch - Mean Velocity', 'Snatch - Peak Velocity'});
disp(correlation);
%% kmeans clustering

% Power Clean - Mean Velocity
subplot(2,2,1)
[idx,C] = kmeans(powerclean_meandata(:,6:7),3); % Use the 6th and 7th column, find 3 cluster
% figure
gscatter(powerclean_meandata(:,6),powerclean_meandata(:,7),idx,'gbr'); lsline;
hold on;
plot(C(:,1),C(:,2),'kx','LineWidth',2); title('Power Clean - Mean Velocity (All %1RM)'); xlabel('Mean of FLEX and VICON'); ylabel('Difference between FLEX and VICON');
legend('Cluster 1','Cluster 2', 'Cluster 3')
% figure
% silhouette(powerclean_meandata(:,6:7),idx)
cluster_correlation = zeros(3,2);
for i=1:3
    data2 = [idx powerclean_meandata(:,6:7)];
    data2 = data2(data2(:,1)==i,:); % Trim data to only if column 1 == 1
    [R,P] = corrcoef(data2(:,2),data2(:,3));
    cluster_correlation(i,1) = R(2);
    cluster_correlation(i,2) = P(2);
end
disp('Power Clean - Mean Velocity Cluster Correlation')
cluster_correlation = array2table(cluster_correlation, 'VariableNames', {'r value', 'p value'}); disp(cluster_correlation);

% Power Clean - Peak Velocity
subplot(2,2,2)
[idx,C] = kmeans(powerclean_peakdata(:,6:7),3);
% figure
gscatter(powerclean_peakdata(:,6),powerclean_peakdata(:,7),idx,'gbr'); lsline;
hold on
plot(C(:,1),C(:,2),'kx','LineWidth',2); title('Power Clean - Peak Velocity (All %1RM)'); xlabel('Mean of FLEX and VICON'); ylabel('Difference between FLEX and VICON');
legend('Cluster 1','Cluster 2', 'Cluster 3')
cluster_correlation = zeros(3,2);
for i=1:3
    data2 = [idx powerclean_peakdata(:,6:7)];
    data2 = data2(data2(:,1)==i,:); % Trim data to only if column 1 == 1
    [R,P] = corrcoef(data2(:,2),data2(:,3));
    cluster_correlation(i,1) = R(2);
    cluster_correlation(i,2) = P(2);
end
disp('Power Clean - Peak Velocity Cluster Correlation')
cluster_correlation = array2table(cluster_correlation, 'VariableNames', {'r value', 'p value'}); disp(cluster_correlation);

% Snatch - Mean Velocity
subplot(2,2,3)
[idx,C] = kmeans(snatch_meandata(:,6:7),3);
% figure
gscatter(snatch_meandata(:,6),snatch_meandata(:,7),idx,'gbr'); lsline;
hold on
plot(C(:,1),C(:,2),'kx','LineWidth',2); title('Snatch - Mean Velocity (All %1RM)'); xlabel('Mean of FLEX and VICON'); ylabel('Difference between FLEX and VICON');
legend('Cluster 1','Cluster 2', 'Cluster 3')
cluster_correlation = zeros(3,2);
for i=1:3
    data2 = [idx snatch_meandata(:,6:7)];
    data2 = data2(data2(:,1)==i,:); % Trim data to only if column 1 == 1
    [R,P] = corrcoef(data2(:,2),data2(:,3));
    cluster_correlation(i,1) = R(2);
    cluster_correlation(i,2) = P(2);
end
disp('Snatch - Mean Velocity Cluster Correlation')
cluster_correlation = array2table(cluster_correlation, 'VariableNames', {'r value', 'p value'}); disp(cluster_correlation);

% Snatch - Peak Velocity 
subplot(2,2,4)
[idx,C] = kmeans(snatch_peakdata(:,6:7),3);
% figure
gscatter(snatch_peakdata(:,6),snatch_peakdata(:,7),idx,'gbr'); lsline;
hold on
plot(C(:,1),C(:,2),'kx','LineWidth',2); title('Snatch - Peak Velocity (All %1RM)'); xlabel('Mean of FLEX and VICON'); ylabel('Difference between FLEX and VICON');
legend('Cluster 1','Cluster 2', 'Cluster 3')
cluster_correlation = zeros(3,2);
for i=1:3
    data2 = [idx snatch_peakdata(:,6:7)];
    data2 = data2(data2(:,1)==i,:); % Trim data to only if column 1 == 1
    [R,P] = corrcoef(data2(:,2),data2(:,3));
    cluster_correlation(i,1) = R(2);
    cluster_correlation(i,2) = P(2);
end
disp('Snatch - Peak Velocity Cluster Correlation')
cluster_correlation = array2table(cluster_correlation, 'VariableNames', {'r value', 'p value'}); disp(cluster_correlation);
%% Calculate CV = SD/Mean

% Column 4 = FLEX
cv_flex = std(data(:,4))/mean(data(:,4))*100;
% Column 5 = VICON
cv_vicon = std(data(:,5))/mean(data(:,5))*100;

%%
tscore = tinv(.95,size(data,2));
standarderror = std(data(:,7))/sqrt(size(data,2));
upperconfidence = meandifferences + tscore * standarderror;
lowerconfidence = meandifferences - tscore * standarderror;
%% 
toc