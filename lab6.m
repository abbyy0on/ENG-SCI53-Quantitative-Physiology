%%
% Author: Abby Yoon
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

%% SPO2 and HR

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