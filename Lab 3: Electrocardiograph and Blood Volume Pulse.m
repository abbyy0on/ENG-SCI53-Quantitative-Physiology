% load all subjects
subject1 = load('Lab3_Ex2_Sub1.mat');
subject2 = load('Lab3_Ex2_Sub2.mat');
subject3 = load('Lab3_Ex2_Sub3.mat');

%% sitting data 

% extract sitting data for all subjects
subject1_sitting = subject1.data(subject1.datastart(1, 1) : subject1.dataend(1,1));
subject2_sitting = subject2.data(subject2.datastart(1, 1) : subject2.dataend(1,1));
subject3_sitting = subject3.data(subject3.datastart(1, 1) : subject3.dataend(1,1));

% subject1 peaks
e20_1 = subject1_sitting / max(subject1_sitting);
de20_1 = diff(e20_1); % derivative of this vector
threshold = 0.8; % peak treshold
ind_1 = find((e20_1(2:end-1)>threshold)&(de20_1(1:end-1)>0)&(de20_1(2:end)<0))+1; % the indices of all peaks in the first 20 second segment

% subject2 peaks
e20_2 = subject2_sitting / max(subject2_sitting);
de20_2 = diff(e20_2); % derivative of this vector
threshold = 0.8; % peak treshold
ind_2 = find((e20_2(2:end-1)>threshold)&(de20_2(1:end-1)>0)&(de20_2(2:end)<0))+1;

% subject3 peaks
e20_3 = subject3_sitting / max(subject3_sitting);
de20_3 = diff(e20_3); % derivative of this vector
threshold = 0.8; % peak treshold
ind_3 = find((e20_3(2:end-1)>threshold)&(de20_3(1:end-1)>0)&(de20_3(2:end)<0))+1;

% find average of peaks / R  waves 
meanSittingSub1 = mean(subject1_sitting(ind_1));
meanSittingSub2 = mean(subject2_sitting(ind_2));
meanSittingSub3 = mean(subject3_sitting(ind_3));

%% standing data

% extract standing data for all subjects
subject1_standing = subject1.data(subject1.datastart(1, 2) : subject1.dataend(1,2));
subject2_standing = subject2.data(subject2.datastart(1, 2) : subject2.dataend(1,2));
subject3_standing = subject3.data(subject3.datastart(1, 2) : subject3.dataend(1,2));

% subject1 peaks
f20_1 = subject1_standing / max(subject1_standing);
df20_1 = diff(f20_1); % derivative of this vector
threshold = 0.8; % peak treshold
ind_11 = find((f20_1(2:end-1)>threshold)&(df20_1(1:end-1)>0)&(df20_1(2:end)<0))+1; % the indices of all peaks in the first 20 second segment

% subject2 peaks
f20_2 = subject2_standing / max(subject2_standing);
df20_2 = diff(f20_2); % derivative of this vector
threshold = 0.8; % peak treshold
ind_12 = find((f20_2(2:end-1)>threshold)&(df20_2(1:end-1)>0)&(df20_2(2:end)<0))+1; % the indices of all peaks in the first 20 second segment

% subject3 peaks
f20_3 = subject3_standing / max(subject3_standing);
df20_3 = diff(f20_3); % derivative of this vector
threshold = 0.8; % peak treshold
ind_13 = find((f20_3(2:end-1)>threshold)&(df20_3(1:end-1)>0)&(df20_3(2:end)<0))+1; % the indices of all peaks in the first 20 second segment

% find average of peaks / R  waves 

meanStandingSub1 = mean(subject1_standing(ind_11));
meanStandingSub2 = mean(subject2_standing(ind_12));
meanStandingSub3 = mean(subject3_standing(ind_13));

%% post exercise data
% extract exercise data for all subjects 
% NOTE: In Dr. Moyer's data, only Subject 1 had any data on post-exercise
% activity
subject1_exercise = subject1.data(subject1.datastart(1, 3) : subject1.dataend(1,3));
 
% subject1 peaks
g20_1 = subject1_exercise / max(subject1_exercise);
dg20_1 = diff(g20_1); % derivative of this vector
threshold = 0.8; % peak treshold
ind_111 = find((g20_1(2:end-1)>threshold)&(dg20_1(1:end-1)>0)&(dg20_1(2:end)<0))+1; % the indices of all peaks in the first 20 second segment

% generate summary for all data
meanExerciseSub1 = mean(subject1_exercise(ind_111))
stdExerciseSub1 = std(subject1_exercise(ind_111))

%% plot

% plot 

std_dev = [stdSittingSub1 stdSittingSub2 stdSittingSub3;
           stdStandingSub1 stdStandingSub2 stdStandingSub3];

figure 
y = [meanSittingSub1 meanSittingSub2 meanSittingSub3;
   meanStandingSub1 meanStandingSub2 meanStandingSub3 ];
bar1 = bar(y);
% set the x-axis labels
x = ["Sitting" "Standing"];
xticks(1:numel(x));  % set the x-ticks at positions 1 and 2
xticklabels(x);     % set the x-tick labels to the categories in 'x'
title("Average Amplitude for R Wave Across Team Members")
ylabel("Voltage (V)")
xlabel("Position")
legend("Subject 1", "Subject 2", "Subject 3")

hold on 

[groups, bars] = size(y);
x = nan(bars, groups);
for i = 1:bars
    x(i,:) = bar1(i).XEndPoints;
end

errorbar(x', y, std_dev, 'k', 'linestyle', 'none', 'HandleVisibility', 'off', "LineWidth", 1.5)
hold off

 
