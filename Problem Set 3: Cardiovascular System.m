%%
% File Name: pset.3
% Author: Abby Yoon
% Created: Oct 2, 2023
% Description: P-set 3 Answers

clear % Clears the workspace variables and Command Window
close all % Close all figure windows
clc % Clears your Command Window

%% 6  

ecg = load('ECGdata.mat');

% no smoothing
figure
subplot(2,2,1)
plot(ecg.ptime, ecg.pdata)
title("no smoothing");
ylabel("ECG data");
xlabel("time");

% plot the same ECG, smoothed, with a span of 10 points

subplot(2,2,2)
ecg.pdata10 = smooth(ecg.pdata,10);
plot(ecg.ptime, ecg.pdata10);
title("smoothed span = 10");
ylabel("ECG data");
xlabel("time");

% plot with a smoothed span of 100 points
subplot(2,2,3)
ecg.pdata100 = smooth(ecg.pdata,100);
plot(ecg.ptime, ecg.pdata100);
title("smoothed span = 100");
ylabel("ECG data");
xlabel("time");


% plot with a smoothed span of 1000 points
subplot(2,2,4)
ecg.pdata1000 = smooth(ecg.pdata,1000);
plot(ecg.ptime, ecg.pdata1000);
title("smoothed span = 1000");
ylabel("ECG data");
xlabel("time");

% The advantages of smoothing is that it decreases noise shown in the
% graph and improves the readability of the figure. The disadvantages are
% that it can remove outliers or stark trends and numbers that are
% typically removed in order to reduce 'noise.' If the graph is smoothed so
% much, it can distort the data and take away from its true values.

% new figure w/ peaks and trouphs 
ecg.pdata30 = smooth(ecg.pdata,30); % Dr. Moyer in Office Hours said to change the smoothing span, 100 is too large, 10 is too small; QRS waves should be clearly seen 
dy = diff(ecg.pdata30);
signDY  =  sign(dy);  
diffsignd = diff(signDY);

% peaks
thresholdMax = 0.65*max(ecg.pdata30); % set threshold so that matLab doesn't find intermediate peaks/troughs!
peaks = find(diffsignd == -2); % an array w/ indices where signDY switches from 1 to -1 (peaks)
peaksData = ecg.pdata30(peaks); 
relevantPeaksData = find(peaksData >= thresholdMax);
indicesRelevantPeaks = peaks(relevantPeaksData);

% Troughs
thresholdMin = 0.65*min(ecg.pdata30);
troughs = find(diffsignd == 2);
troughsData = ecg.pdata30(troughs);
relevantTroughsData = find(troughsData <= thresholdMin);
indicesRelevantTroughs = troughs(relevantTroughsData);

figure 
hold on 
plot(ecg.ptime, ecg.pdata30)
plot(ecg.ptime(indicesRelevantPeaks), ecg.pdata30(indicesRelevantPeaks), "ys", 'linewidth', 2.5);
plot(ecg.ptime(indicesRelevantTroughs), ecg.pdata30(indicesRelevantTroughs), "ys",'linewidth', 2.5);
title("Peaks and Troughs of Best Smoothed ECG Trace");
ylabel("ECG Data");
xlabel("Time (t)");
legend("ECG Data", "Peaks and Troughs");
hold off
