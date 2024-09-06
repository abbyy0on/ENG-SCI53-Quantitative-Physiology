%%
% Author: Abby Yoon

%% Figure 1
abby_datastart = abby_data.datastart;
adam_datastart = adam_data.datastart;
alice_datastart = alice_data.datastart;

% datastart(channel #, block # (2nd column) of comment #)+ comtickpos(row #)
ay_comtickpos = abby_data.com(:,3); % look at 3rd column
aw_comtickpos = adam_data.com(:,3);
ac_comtickpos = alice_data.com(:,3);
samplerate = 100; %ticks/second
ay_time = (1/(abby_data.samplerate(1)):1/(abby_data.samplerate(1)):60);
aw_time = (1/(adam_data.samplerate(1)):1/(adam_data.samplerate(1)):60);
ac_time = (1/(alice_data.samplerate(1)):1/(alice_data.samplerate(1)):60);
ay_shifted_time = ((1/(abby_data.samplerate(1))):1/
(abby_data.samplerate(1)):60+200/samplerate);

% 1st comment = baseline
abby_normal = abby_data.data(abby_datastart(1,1)+ay_comtickpos(1):
(abby_datastart(1,1)+ay_comtickpos(1)+ 60*samplerate(1))-1);
adam_normal = adam_data.data(adam_datastart(1,1)+aw_comtickpos(1):
(adam_datastart(1,1)+aw_comtickpos(1)+ 60*samplerate(1))-1);
alice_normal = alice_data.data(alice_datastart(1,1)+ac_comtickpos(1):
(alice_datastart(1,1)+ac_comtickpos(1)+ 60*samplerate(1))-1);

% 2nd comment = inhale
% moved 200 ticks to the left for abby's data for better view
abby_inhale = abby_data.data(abby_datastart(1,1)+ay_comtickpos(2)-200:
(abby_datastart(1,1)+ay_comtickpos(2)+60*samplerate(1))-1);
adam_inhale = adam_data.data(adam_datastart(1,1)+aw_comtickpos(2):
(adam_datastart(1,1)+aw_comtickpos(2)+60*samplerate(1))-1);
alice_inhale = alice_data.data(alice_datastart(1,1)+ac_comtickpos(2):
(alice_datastart(1,1)+ac_comtickpos(2)+60*samplerate(1))-1);

% 4th comment = exhale
abby_exhale = abby_data.data(abby_datastart(1,1)+ay_comtickpos(4)-200:
(abby_datastart(1,1)+ay_comtickpos(4)+60*samplerate(1))-1);
adam_exhale = adam_data.data(adam_datastart(1,1)+aw_comtickpos(4):
(adam_datastart(1,1)+aw_comtickpos(4)+60*samplerate(1))-1);
alice_exhale = alice_data.data(alice_datastart(1,1)+ac_comtickpos(4):
(alice_datastart(1,1)+ac_comtickpos(4)+60*samplerate(1))-1);

%% Plot fig 1
subplot(3, 3, 1);
plot(ay_time, abby_normal);
title('Abby - Normal Breathing');
xlabel('Time (s)');
ylabel('Voltage (mV)');
subplot(3, 3, 4);
plot(ay_shifted_time, abby_inhale);
title('Abby - Inhale and Hold');
xlabel('Time (s)');
ylabel('Voltage (mV)');
subplot(3, 3, 7);
plot(ay_shifted_time, abby_exhale);
title('Abby - Exhale and Hold');
xlabel('Time (s)');
ylabel('Voltage (mV)');

% Adam
subplot(3, 3, 2);
plot(aw_time, adam_normal);
title('Adam - Normal Breathing');
xlabel('Time (s)');
ylabel('Voltage (mV)');
subplot(3, 3, 5);
plot(aw_time, adam_inhale);
title('Adam - Inhale and Hold');
xlabel('Time (s)');
ylabel('Voltage (mV)');
subplot(3, 3, 8);
plot(aw_time, adam_exhale);
title('Adam - Exhale and Hold');
xlabel('Time (s)');
ylabel('Voltage (mV)');

% Alice
subplot(3, 3, 3);
plot(ac_time, alice_normal);
title('Alice - Normal Breathing');
xlabel('Time (s)');
ylabel('Voltage (mV)');
subplot(3, 3, 6);
plot(ac_time, alice_inhale);
title('Alice - Inhale and Hold');
xlabel('Time (s)');
ylabel('Voltage (mV)');
subplot(3, 3, 9);
plot(ac_time, alice_exhale);
title('Alice - Exhale and Hold');
xlabel('Time (s)');
ylabel('Voltage (mV)');

%% Table 1
calculate_breathing = @(data, samplerate) length(findpeaks(data))/
(length(data)/samplerate);
ay_br_per_min = calculate_breathing(abby_normal, samplerate);
aw_br_per_min = calculate_breathing(adam_normal, samplerate);
ac_br_per_min = calculate_breathing(alice_normal, samplerate);

abby_inhale_breath = ((abby_datastart(1,1)+ay_comtickpos(3))-
(abby_datastart(1,1)+ay_comtickpos(2)))/samplerate;
adam_inhale_breath = ((adam_datastart(1,1)+aw_comtickpos(3))-
(adam_datastart(1,1)+aw_comtickpos(2)))/samplerate;
alice_inhale_breath = ((alice_datastart(1,1)+ac_comtickpos(3))-
(alice_datastart(1,1)+ac_comtickpos(2)))/samplerate;

%exhalation hold
abby_exhale_hold = ((abby_datastart(1,1)+ay_comtickpos(5))-
(abby_datastart(1,1)+ay_comtickpos(4)))/samplerate;
adam_exhale_hold = ((adam_datastart(1,1)+aw_comtickpos(5))-
(adam_datastart(1,1)+aw_comtickpos(4)))/samplerate;
alice_exhale_hold = ((alice_datastart(1,1)+ac_comtickpos(5))-
(alice_datastart(1,1)+ac_comtickpos(4)))/samplerate;

%% Table 2
abby_hr_normal = abby_data.data(abby_datastart(3,1)+ay_comtickpos(1):
(abby_datastart(3,1)+ay_comtickpos(1)+ 60*samplerate(1))-1);
adam_hr_normal = adam_data.data(adam_datastart(3,1)+aw_comtickpos(1):
(adam_datastart(3,1)+aw_comtickpos(1)+ 60*samplerate(1))-1);
alice_hr_normal = alice_data.data(alice_datastart(3,1)+ac_comtickpos(1):
 (alice_datastart(3,1)+ac_comtickpos(1) +60*samplerate(1))-1);
abby_hr_mean_normal = mean(abby_hr_normal);
adam_hr_mean_normal = mean(adam_hr_normal);
alice_hr_mean_normal = mean(alice_hr_normal);

abby_hr_inhale = abby_data.data(abby_datastart(3,1)+ay_comtickpos(2):
(abby_datastart(3,1)+ay_comtickpos(3)+ samplerate(1))-1);
adam_hr_inhale = adam_data.data(adam_datastart(3,1)+aw_comtickpos(2):
(adam_datastart(3,1)+aw_comtickpos(3)+ samplerate(1))-1);
alice_hr_inhale = alice_data.data(alice_datastart(3,1)+ac_comtickpos(2):
(alice_datastart(3,1)+ac_comtickpos(3)+ samplerate(1))-1);
abby_hr_mean_inhale = mean(abby_hr_inhale);
adam_hr_mean_inhale = mean(adam_hr_inhale);
alice_hr_mean_inhale = mean(alice_hr_inhale);

abby_hr_exhale = abby_data.data(abby_datastart(3,1)+ay_comtickpos(4):
(abby_datastart(3,1)+ay_comtickpos(5)+ samplerate(1))-1);
adam_hr_exhale = adam_data.data(adam_datastart(3,1)+aw_comtickpos(4):
(adam_datastart(3,1)+aw_comtickpos(5)+ samplerate(1))-1);
alice_hr_exhale = alice_data.data(alice_datastart(3,1)+ac_comtickpos(4):
(alice_datastart(3,1)+ac_comtickpos(5)+ samplerate(1))-1);
abby_hr_mean_exhale = mean(abby_hr_exhale);
adam_hr_mean_exhale = mean(adam_hr_exhale);
alice_hr_mean_exhale = mean(alice_hr_exhale);

%% Exercise 2  
%% Figure 2
% load the data file for all exercises & subjects
abby_data = load('Lab6_Abby.mat');
alice_data = load('Lab6_Alice.mat');
adam_data = load('Lab6_Adam.mat');

% the sample rate is the same across all blocks, channels, & group members
sr = abby_data.samplerate(1); % indexes/extracts the first element 
dt = 1/sr;

% extract ex2 data for abby
abby_ex2_breath = abby_data.data(abby_data.datastart(1, 2) : abby_data.dataend(1,2));
abby_ex2_pulse = abby_data.data(abby_data.datastart(2, 2) : abby_data.dataend(2,2));
abby_ex2_hr = abby_data.data(abby_data.datastart(3, 2) : abby_data.dataend(3,2));

% extract ex2 data for alice
alice_ex2_breath = alice_data.data(alice_data.datastart(1, 2) : alice_data.dataend(1,2));
alice_ex2_pulse = alice_data.data(alice_data.datastart(2, 2) : alice_data.dataend(2,2));
alice_ex2_hr = alice_data.data(alice_data.datastart(3, 2) : alice_data.dataend(3,2));

% extract ex2 data for alice
adam_ex2_breath = adam_data.data(adam_data.datastart(1, 2) : adam_data.dataend(1,2));
adam_ex2_pulse = adam_data.data(adam_data.datastart(2, 2) : adam_data.dataend(2,2));
adam_ex2_hr = adam_data.data(adam_data.datastart(3, 2) : adam_data.dataend(3,2));

 %%  SPO2 and HR
% read SPO2 data
spo2 = 'ES 53 Lab 6 - SPO2 Data.xlsx';

% ABBY
% make a time vector
time_abby = (1:length(abby_ex2_breath)) * dt;  

% read appropriate data in excel
abby_data_table = readtable(spo2, 'Sheet', 'Abby');  
% extract two columns
abby_spo2_ex2 = abby_data_table.SPO2_1.'; 
abby_time_ex2 = abby_data_table.Time_1.'

% find heartrate
[peaks, locs] = findpeaks(abby_ex2_pulse,'MinPeakHeight',.01,'MinPeakProminence',.0001); % threshold was set to .01 to obtain highest peaks
locs_abby = locs; % rename 'locs' and 'peaks' to be specific to each group member + stay organized
peaks_abby = peaks;

figure % verifying  correct peaks  
plot(abby_ex2_pulse)
hold on
plot(locs_abby,peaks_abby,"x")
d_abby=diff(locs_abby);
heartrate_abby = 60.*d_abby.*dt
figure   %  verifying  correct HR  data
plot(heartrate_abby)


%  ALICE
% make a time vector
time_alice = (1:length(alice_ex2_breath)) * dt;  

% read appropriate data in excel
alice_data_table = readtable(spo2, 'Sheet', 'Alice');  

% extract two columns
alice_spo2_ex2 = alice_data_table.SPO2_1.'; 
alice_time_ex2 = alice_data_table.Time_1.'

% find heartrate
[peaks, locs] = findpeaks(alice_ex2_pulse,'MinPeakHeight',.014,'MinPeakProminence',.0001);
locs_alice = locs; % rename 'locs' and 'peaks' to be specific to each group member + stay organized
peaks_alice = peaks;

figure % verifying  correct peaks  
plot(alice_ex2_pulse)
hold on
plot(locs_alice,peaks_alice,"x")
d_alice=diff(locs_alice);
heartrate_alice = 60.*d_alice.*dt
figure   %  verifying  correct HR  data
plot(heartrate_alice)


%  ADAM
% make a time vector
time_adam = (1:length(adam_ex2_breath)) * dt;  

% read appropriate data in excel
adam_data_table = readtable(spo2, 'Sheet', 'Adam');  

% extract two columns
adam_spo2_ex2 = adam_data_table.SPO2_1.'; 
adam_time_ex2 = adam_data_table.Time_1.'

% find heartrate
[peaks, locs] = findpeaks(adam_ex2_pulse,'MinPeakHeight',.006,'MinPeakProminence',.0001);
locs_adam = locs; % rename 'locs' and 'peaks' to be specific to each group member + stay organized
peaks_adam = peaks;

figure % verifying  correct peaks  
plot(adam_ex2_pulse)
hold on
plot(locs_adam,peaks_adam,"x")
d_adam=diff(locs_adam);
heartrate_adam = 60.*d_adam.*dt
figure   %  verifying  correct HR  data
plot(heartrate_adam)

%% Make Subplots

% ABBY
% make sure length of  vector  aligns  for manual HR data
time_pulse_abby = locs_abby(1:end-1)*dt

% make subplots
figure
subplot(3,1,1)

% Plot Breath Signal
yyaxis left
plot(time_abby, abby_ex2_breath, 'b');
xlabel("Time (s)");
ylabel("voltage (mV)");

% Plot SPO2 data
hold on
plot(abby_time_ex2, abby_spo2_ex2, '-m');
ylabel("SPO2 (%)");

yyaxis right
% Plot Heart Rate data on the right y-axis
plot(time_pulse_abby, heartrate_abby, 'r'); 
ylabel("Heart Rate (bpm)");

% Add legend
legend("Breathing Rate", "SPO2", "Heart Rate");
legend('FontSize', 7)


% ADAM
% make sure length of  vector  aligns  for manual HR data
time_pulse_adam = locs_adam(1:end-1)*dt

% make subplots
subplot(3,1,2)

% Plot Breath Signal
yyaxis left
plot(time_adam, adam_ex2_breath, 'b');
xlabel("Time (s)");
ylabel("voltage (mV)");

% Plot SPO2 data
hold on
plot(adam_time_ex2, adam_spo2_ex2, '-m');
ylabel("SPO2 (%)");

yyaxis right
% Plot Heart Rate data on the right y-axis
plot(time_pulse_adam, heartrate_adam, 'r'); 
ylabel("Heart Rate (bpm)");

% Add legend
legend("Breathing Rate", "SPO2", "Heart Rate");
legend('FontSize', 7)


% ALICE
% make sure length of  vector  aligns  for manual HR data
time_pulse_alice = locs_alice(1:end-1)*dt

% make subplots
subplot(3,1,3)

% Plot Breath Signal
yyaxis left
plot(time_alice, alice_ex2_breath, 'b');
xlabel("Time (s)");
ylabel("voltage (mV)");

% Plot SPO2 data
hold on
plot(alice_time_ex2, alice_spo2_ex2, '-m');
ylabel("SPO2 (%)");

yyaxis right
% Plot Heart Rate data on the right y-axis
plot(time_pulse_alice, heartrate_alice, 'r'); 
ylabel("Heart Rate (bpm)");

% Add legend
legend("Breathing Rate", "SPO2", "Heart Rate");
legend('FontSize', 7)

%% Table 3 - average breathing rate during hyperventilation

abby_hypervent_data = abby_data.data(abby_data.datastart(1,2) + abby_data.com(7,3) : abby_data.datastart(1,2) + abby_data.com(8,3)); % 'hyperventilate' was our 7th comment; 'breathe' was our 8th comment
alice_hypervent_data = alice_data.data(alice_data.datastart(1,2) + alice_data.com(7,3) : alice_data.datastart(1,2) + alice_data.com(8,3));
adam_hypervent_data = adam_data.data(adam_data.datastart(1,2) + adam_data.com(7,3) : adam_data.datastart(1,2) + adam_data.com(8,3));

% smooth data
abby_hypervent = smooth(abby_hypervent_data, 40)
adam_hypervent = smooth(adam_hypervent_data, 40)
alice_hypervent = smooth(alice_hypervent_data, 40)

[abby_hventpeaks,abby_hventlocs] = findpeaks(abby_hypervent,'MinPeakDistance', 2);  
[adam_hventpeaks,adam_hventlocs] = findpeaks(adam_hypervent,'MinPeakHeight', 50,'MinPeakDistance', 2);  
[alice_hventpeaks,alice_hventlocs] = findpeaks(alice_hypervent,'MinPeakHeight', 70,'MinPeakDistance', 2);  

% verify that we have correct peaks
time_abby = (1:length(abby_hypervent)) * dt; 
figure
plot(time_abby, abby_hypervent)
hold on
plot(time_abby(abby_hventlocs), abby_hventpeaks, 'ok')

% verify that we have correct peaks
time_adam = (1:length(adam_hypervent)) * dt; 
figure
plot(time_adam, adam_hypervent)
hold on
plot(time_adam(adam_hventlocs), adam_hventpeaks, 'ok')

%  verify that we have correct peaks
time_alice = (1:length(alice_hypervent)) * dt; 
figure
plot(time_alice, alice_hypervent)
hold on
plot(time_alice(alice_hventlocs), alice_hventpeaks, 'ok')

% find breathing rate calculations
abby_diff = diff(abby_hventlocs);
alice_diff = diff(alice_hventlocs);
adam_diff = diff(adam_hventlocs);
abby_breathing = 60.*abby_diff.*dt;
alice_breathing = 60.*alice_diff.*dt;
adam_breathing = 60.*adam_diff.*dt;

abby_brate = mean(abby_breathing)
alice_brate = mean(alice_breathing)
adam_brate = mean(adam_breathing)

% find period of hyperventilation
abby_hvent_period = length(time_abby)*dt
alice_hvent_period = length(time_alice)*dt
adam_hvent_period = length(time_adam)*dt

% find initial + final SPO2 for abby
time_before_hvent_abby = (abby_data.com(7,3)-1)*dt % seconds passed in ex2 until 'hyperventilation' comment
time_after_hvent_abby = (abby_data.com(8,3)+1)*dt % seconds passed in ex2 until after 'breathe' comment
start_ex2_abby = (abby_data.com(6,3))*dt % time at which ex2 began
time_prehvent_abby_spo2 = (time_before_hvent_abby - start_ex2_abby) % calculate how many seconds passed from the beginning of ex 2 until hyperventilation
time_prehvent_abby_spo2 = round(time_prehvent_abby_spo2 / 10) * 10 % round to nearest multiple of 10
time_posthvent_abby_spo2 = (time_after_hvent_abby - start_ex2_abby) % calculate how many seconds passed from the beginning of ex 2 until hyperventilation ended
time_posthvent_abby_spo2 = round(time_posthvent_abby_spo2 / 10) * 10 % round to nearest multiple of 10
abby_spo2_ex2_combined = [abby_time_ex2', abby_spo2_ex2']; % combine time & SPO2 data into one vector
spo2_initial_abby = abby_spo2_ex2_combined(abby_spo2_ex2_combined(:,1) == time_prehvent_abby_spo2, 2) % find the spo2 that corresponds w/ time of interest
spo2_final_abby = abby_spo2_ex2_combined(abby_spo2_ex2_combined(:,1) == time_posthvent_abby_spo2, 2) % find the spo2 that corresponds w/ time of interest

% find initial + final SPO2 for adam
time_before_hvent_adam = (adam_data.com(7,3)-1)*dt % seconds passed in ex2 until 'hyperventilation' comment
time_after_hvent_adam = (adam_data.com(8,3)+1)*dt % seconds passed in ex2 until after 'breathe' comment
start_ex2_adam = (adam_data.com(6,3))*dt % time at which ex2 began
time_prehvent_adam_spo2 = (time_before_hvent_adam - start_ex2_adam) % calculate how many seconds passed from the beginning of ex 2 until hyperventilation
time_prehvent_adam_spo2 = round(time_prehvent_adam_spo2 / 10) * 10 % round to nearest multiple of 10
time_posthvent_adam_spo2 = (time_after_hvent_adam - start_ex2_adam) % calculate how many seconds passed from the beginning of ex 2 until hyperventilation ended
time_posthvent_adam_spo2 = round(time_posthvent_adam_spo2 / 10) * 10 % round to nearest multiple of 10
adam_spo2_ex2_combined = [adam_time_ex2', adam_spo2_ex2']; % combine time & SPO2 data into one vector
spo2_initial_adam = adam_spo2_ex2_combined(adam_spo2_ex2_combined(:,1) == time_prehvent_adam_spo2, 2) % find the spo2 that corresponds w/ time of interest
spo2_final_adam = adam_spo2_ex2_combined(adam_spo2_ex2_combined(:,1) == time_posthvent_adam_spo2, 2) % find the spo2 that corresponds w/ time of interest

% find initial + final SPO2 for alice
time_before_hvent_alice = (alice_data.com(7,3)-1)*dt % seconds passed in ex2 until 'hyperventilation' comment
time_after_hvent_alice = (alice_data.com(8,3)+1)*dt % seconds passed in ex2 until after 'breathe' comment
start_ex2_alice = (alice_data.com(6,3))*dt % time at which ex2 began
time_prehvent_alice_spo2 = (time_before_hvent_alice - start_ex2_alice) % calculate how many seconds passed from the beginning of ex 2 until hyperventilation
time_prehvent_alice_spo2 = round(time_prehvent_alice_spo2 / 10) * 10 % round to nearest multiple of 10
time_posthvent_alice_spo2 = (time_after_hvent_alice - start_ex2_alice) % calculate how many seconds passed from the beginning of ex 2 until hyperventilation ended
time_posthvent_alice_spo2 = round(time_posthvent_alice_spo2 / 10) * 10 % round to nearest multiple of 10
alice_spo2_ex2_combined = [alice_time_ex2', alice_spo2_ex2']; % combine time & SPO2 data into one vector
spo2_initial_alice = alice_spo2_ex2_combined(alice_spo2_ex2_combined(:,1) == time_prehvent_alice_spo2, 2) % find the spo2 that corresponds w/ time of interest
spo2_final_alice = alice_spo2_ex2_combined(alice_spo2_ex2_combined(:,1) == time_posthvent_alice_spo2, 2) % find the spo2 that corresponds w/ time of interest

% find initial + final HR for abby
initial_hr_abby = abby_ex2_hr(time_before_hvent_abby/dt)
final_hr_abby = abby_ex2_hr(time_after_hvent_abby/dt)

% find initial + final HR for adam
initial_hr_adam = adam_ex2_hr(time_before_hvent_adam/dt)
final_hr_adam = adam_ex2_hr(time_after_hvent_adam/dt)

% find initial + final HR for alice
initial_hr_alice = alice_ex2_hr(time_before_hvent_alice/dt)
final_hr_alice = alice_ex2_hr(time_after_hvent_alice/dt)

%% Exercise 3
data_abby = load("Lab6_Abby.mat");
data_adam = load("Lab6_Adam.mat");
data_alice = load("Lab6_Alice.mat");

spo2_abby = readmatrix("ES 53 Lab 6 - SPO2 Data.xlsx", "Sheet", "Abby");
spo2_adam = readmatrix("ES 53 Lab 6 - SPO2 Data.xlsx", "Sheet", "Adam");
spo2_alice = readmatrix("ES 53 Lab 6 - SPO2 Data.xlsx", "Sheet", "Alice");

%% Figure 3 and Table 4
data = [data_abby, data_adam, data_alice];
spo2_data = {spo2_abby, spo2_adam, spo2_alice};
blc = 3;
chn = 1;
cn = 0;
from = 0;
to = 240;
hr_thresholds = [0.015 0.005 0.01];
br_thresholds = [50 18 56];
table_4 = zeros(6,3);
for i=1:length(data)
 [t, breath] = Extract(data(i), blc, chn, cn, from, to);
 % Extract pulse to find heart rate
 [~, pulse] = Extract(data(i), blc, 2, cn, from, to);
 [t_hr, hr] = Heartrate(t, pulse, hr_thresholds(i));
 %[~, lab_hr] = Extract(data(i), blc, 3, cn, from, to);
 spo2 = spo2_data{i};
spo2 = spo2(1:25,7); % Take first 25 timepts (240s)
 t_spo2 = (0:length(spo2)-1) .* 10;
 figure(3);
 subplot(3,1,i);
 plot(t, breath, 'b');
 xlim([0 240])
 ylabel("Voltage (mV)");
 xlabel("Time (s)");
 addaxis(t_hr, hr, 'r');
 addaxislabel(2, "Heart Rate (bpm)");
 addaxis(t_spo2, spo2);
 addaxislabel(3, "SPO_2 (%)")
 legend("Breathing Rate", "Heart Rate", "SPO_2");
 % Make rebreathing table
 rbstart = data(i).com(13,3); % Rebreathing start index
 rbend = data(i).com(14,3);
 % Average rebreathing rate
 [t_br, br] = Breathrate(t, breath, br_thresholds(i));
 table_4(1,i) = mean(br(find(t_br > t(rbstart) & t_br < t(rbend))));
 table_4(2,i) = t(rbend) - t(rbstart); % Rebreathing period
 table_4(3,i) = spo2(floor(t(rbstart)/10)+1); % Initial SPO2 percentage
 table_4(4,i) = spo2(floor(t(rbend)/10)+1); % Final SPO2 percentage
 table_4(5,i) = hr(find(t_hr > t(rbstart) - 1, 1)); % Initial heart rate
 table_4(6,i) = hr(find(t_hr > t(rbend) - 1, 1)); % Final heart rate
end
disp("Table 4");
disp(table_4);

%% Extraction function
% Extract portion of data starting at a specified comment or time through
% 'end' seconds later or through end if 'to' is 0
function [t,snip] = Extract(b,blc,chn,cn,from,to)
 %starting index point of snip of data that I want
 if cn == 0
 comint = 0;
 else
 comint = b.com(cn,3);
 end
 bfrom = b.datastart(chn,blc)+comint+from*b.samplerate(1);
 if to == -1 %specify interval = -1 if you want entire block
 bto = b.dataend(chn,blc);
 elseif to == 0 %specify interval = 0 if you want interval to next comment
 bto = b.datastart(chn,blc)+b.com(cn+1,3);
else
 % Ending index of point of snip of data that I want
 bto = b.datastart(chn,blc)+comint+to*b.samplerate(1)-1;
 end
 snip = b.data(bfrom:bto); %select data between "from" and "to"
 t = (0:length(snip)-1)/b.samplerate(1);
end

%% Heart rate function
function [t_hr, hr] = Heartrate(t, data, threshold)
 samplerate = 1/(t(2)-t(1));
 [pks, locs] = findpeaks(data, samplerate, 'MinPeakDistance',0.5, ...
 'MinPeakHeight',threshold);
 % Plot to check thresholds
 %figure;
 %plot(t, data);
 %hold on;
 %scatter(locs, pks);
 %hold off;
 locs = smooth(locs, 15); % Adjust as necessary
 int = diff(locs);
 hr = 1./int .* 60; % bpm
 locs(end) = [];
 t_hr = locs;
end

%% Breathing rate function

function [t_br, br] = Breathrate(t, data, threshold)
 samplerate = 1/(t(2)-t(1));
 [pks, locs] = findpeaks(data, samplerate, 'MinPeakDistance',2, ...
 'MinPeakHeight',threshold);
 % Plot to check thresholds
 %figure;
 %plot(t, data);
 %hold on;
 %scatter(locs, pks);
 %hold off;
 int = diff(locs);
 br = 1./int .* 60;
 locs(end) = [];
 t_br = locs;
end

