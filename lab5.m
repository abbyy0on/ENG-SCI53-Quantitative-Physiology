%% Figure 1 
% load the data file for all exercises & subjects
all_data = load('AAA_ex234.mat');

% extract ex2 FLOW data for all subjects
alice_ex2_flow = all_data.data(all_data.datastart(1, 3) : all_data.dataend(1,3));
abby_ex2_flow = all_data.data(all_data.datastart(1, 10) : all_data.dataend(1,10));
adam_ex2_flow = all_data.data(all_data.datastart(1, 14) : all_data.dataend(1,14));

sr = all_data.samplerate(1); % indexes / extracts the first element 
dt = 1/sr;

% make a time vector
time = (1:length(abby_ex2_flow)); % note: used the array with the least values to avoid errors; abby had the least amount time 

% standardize the lengths of the vectors
alice_ex2_flow = (alice_ex2_flow(1:length(time))); 
abby_ex2_flow = (abby_ex2_flow(1:length(time)));
adam_ex2_flow = (adam_ex2_flow(1:length(time)));

% make subplots
figure
subplot(3,2,1)
plot(time*dt, alice_ex2_flow);
ylim([-2 2.5]); % standardize y axis range for all graphs
ylabel("Flow (L/s)");
xlabel("Time (s)");

subplot(3,2,3)
plot(time*dt, abby_ex2_flow);
ylim([-2 2.5]);
ylabel("Flow (L/s)");
xlabel("Time (s)");

subplot(3,2,5)
plot(time*dt, adam_ex2_flow);
ylim([-2 2.5]);
ylabel("Flow (L/s)");
xlabel("Time (s)");

% extract ex2 VOLUME data for all subjects
alice_ex2_volume = all_data.data(all_data.datastart(2, 3) : all_data.dataend(2,3));
abby_ex2_volume = all_data.data(all_data.datastart(2, 10) : all_data.dataend(2,10));
adam_ex2_volume = all_data.data(all_data.datastart(2, 14) : all_data.dataend(2,14));

sr = all_data.samplerate(1); % indexes / extracts the first element 
dt = 1/sr;

% make a time vector
time = (1:length(abby_ex2_volume)); % note: used the array with the least values to avoid errors; abby had the least amount time 

% standardize the lengths of the vectors
alice_ex2_volume = (alice_ex2_volume(1:length(time))); 
abby_ex2_volume = (abby_ex2_volume(1:length(time)));
adam_ex2_volume = (adam_ex2_volume(1:length(time)));

% perform numerical integration
volume_alice = cumtrapz(time*dt, alice_ex2_flow);  
volume_abby = cumtrapz(time*dt, abby_ex2_flow);  
volume_adam = cumtrapz(time*dt, adam_ex2_flow);  

% make subplots
subplot(3,2,2)
hold on 
plot(time*dt, alice_ex2_volume);
ylabel("Volume (L)");
xlabel("Time (s)");
plot(time*dt, volume_alice)
hold off
ylim([-1 2.5]);
legend("LabChart", "Integrated", 'Location', 'northwest')

subplot(3,2,4)
hold on
plot(time*dt, abby_ex2_volume);
ylabel("Volume (L)");
xlabel("Time (s)");
plot(time*dt, volume_abby)
hold off
ylim([-1 2.5]);
legend("LabChart", "Integrated", 'Location', 'northwest')

subplot(3,2,6)
hold on
plot(time*dt, adam_ex2_volume);
ylabel("Volume (L)");
xlabel("Time (s)");
plot(time*dt, volume_adam)
hold off
ylim([-1 2.5]);
legend("LabChart", "Integrated", 'Location', 'northwest')

%% Table 1 Respiratory Rate
% alice - find peaks in volume data
alice_ex2_volume_smoothed = smooth(alice_ex2_volume,100); % smooth the data so that MatLab doesn't register 'bumps' in the data as peaks
[peaks, locs] = findpeaks(alice_ex2_volume_smoothed);
locs_alice = locs; % rename 'locs' and 'peaks' to be specific to each group member + stay organized
peaks_alice = peaks;

figure % plot the volume vs. time curve and peaks to verify we have correct peaks
plot(time*dt, alice_ex2_volume_smoothed);
hold on
plot(time(locs_alice)*dt, peaks_alice, 'ro'); 
hold off

% number of peaks per minute
peaks_per_minute_alice = length(peaks_alice) / (time(end)*dt)*60

% abby - find peaks
abby_ex2_volume_smoothed = smooth(abby_ex2_volume,100); % smooth the data so that MatLab doesn't register 'bumps' in the data as peaks
[peaks, locs] = findpeaks(abby_ex2_volume_smoothed);
locs_abby = locs;
peaks_abby = peaks;
figure % plot the volume vs. time curve and detected peaks to verify that we have correct peaks
plot(time*dt, abby_ex2_volume_smoothed);
hold on
plot(time(locs_abby)*dt, peaks_abby, 'ro'); 
hold off

% number of peaks per minute
peaks_per_minute_abby = length(peaks_abby) / (time(end)*dt)*60
 
% adam - find peaks
adam_ex2_volume_smoothed = smooth(adam_ex2_volume,100); % smooth the data so that MatLab doesn't register 'bumps' in the data as peaks
[peaks, locs] = findpeaks(adam_ex2_volume_smoothed);
peaks_adam = peaks;
locs_adam = locs;
figure % plot the volume vs. time curve and detected peaks to verify that we have the peaks
plot(time*dt, adam_ex2_volume_smoothed);
hold on
plot(time(locs_adam)*dt, peaks_adam, 'ro'); 
hold off

% number of peaks per minute
peaks_per_minute_adam = length(peaks_adam) / (time(end)*dt)*60

% calculate mean + standard deviation  across all members
mean_peaks_per_minute = mean([peaks_per_minute_adam peaks_per_minute_abby peaks_per_minute_alice])
std_peaks_per_minute = std([peaks_per_minute_adam peaks_per_minute_abby peaks_per_minute_alice])

%% Table 1 Volume of Single Tidal Inspiration  
% alice - find the total volume for all inspirations
positive_flow_indices = alice_ex2_flow > 0;
total_volume_inspiration_alice = trapz(time(positive_flow_indices)*dt, alice_ex2_flow(positive_flow_indices))

% divide this by the number of inspirations that occurred (AKA number of
% peaks) to find the average volume per inspiration
avg_vol_inspiration_alice = total_volume_inspiration_alice / length(peaks_alice)


% abby - find the total volume during all inspirations
positive_flow_indices = abby_ex2_flow > 0;
total_volume_inspiration_abby = trapz(time(positive_flow_indices)*dt,abby_ex2_flow(positive_flow_indices));  
% divide by number of inspirations
avg_vol_inspiration_abby = total_volume_inspiration_abby /  length(peaks_abby)


% adam - find the total volume during all inspirations
positive_flow_indices = adam_ex2_flow > 0;
total_volume_inspiration_adam = trapz(time(positive_flow_indices)*dt,adam_ex2_flow(positive_flow_indices));  
% divide by the number of inspirations 
avg_vol_inspiration_adam = total_volume_inspiration_adam /  length(peaks_adam)

% find mean + std across group members
mean_tidal_volume = mean([avg_vol_inspiration_alice avg_vol_inspiration_abby avg_vol_inspiration_adam])
std_tidal_volume = std([avg_vol_inspiration_alice avg_vol_inspiration_abby avg_vol_inspiration_adam])

%% Table 1 Expired minute volume 
expired_volume_alice = peaks_per_minute_alice * avg_vol_inspiration_alice
expired_volume_abby = peaks_per_minute_abby * avg_vol_inspiration_abby
expired_volume_adam = peaks_per_minute_adam * avg_vol_inspiration_adam

% find mean + std across group members
mean_expired_volume = mean([expired_volume_alice expired_volume_abby expired_volume_adam])
std_expired_volume = std([expired_volume_alice expired_volume_abby expired_volume_adam])

%% Table 1  IRV  
threshold = 1.3; % we want all peaks except the highest peak
medium_indices_adam = find(adam_ex2_volume_smoothed(locs_adam)< threshold); % remove the highest peak
locs_medium_adam = locs_adam(medium_indices_adam); 
figure
plot(time*dt, adam_ex2_volume_smoothed); % plot to verify correct peaks
hold on
plot(time(locs_medium_adam)*dt, adam_ex2_volume_smoothed(locs_medium_adam), "ro");
hold off 

% calculate IRV by finding the average of the medium peaks + subtracting it
% from the max peak
meanMediumPeaks_adam = mean(adam_ex2_volume_smoothed(locs_medium_adam))
irv_adam = max(adam_ex2_volume_smoothed) - meanMediumPeaks_adam


% obtain all peaks except highest peak - abby IRV
medium_indices_abby= find(abby_ex2_volume_smoothed(locs_abby)< threshold); % remove the highest peak
locs_medium_abby = locs_abby(medium_indices_abby); % remove the highest value from peaks_adam
figure
plot(time*dt, abby_ex2_volume_smoothed);
hold on
plot(time(locs_medium_abby)*dt, abby_ex2_volume_smoothed(locs_medium_abby), "ro");
hold off 

% do calculations
meanMediumPeaks_abby = mean(abby_ex2_volume_smoothed(locs_medium_abby))
irv_abby = max(abby_ex2_volume_smoothed) - meanMediumPeaks_abby


% obtain all peaks except highest peak - alice IRV
medium_indices_alice = find(alice_ex2_volume_smoothed(locs_alice)< threshold); % remove the highest peak
locs_medium_alice = locs_alice(medium_indices_alice); 
figure
plot(time*dt, alice_ex2_volume_smoothed);
hold on
plot(time(locs_medium_alice)*dt, alice_ex2_volume_smoothed(locs_medium_alice), "ro");
hold off 

% do calculations
meanMediumPeaks_alice = mean(alice_ex2_volume_smoothed(locs_medium_alice))
irv_alice = max(alice_ex2_volume_smoothed) - meanMediumPeaks_alice

% calculate the mean + std for IRS across all group members
mean_irv = mean([irv_alice irv_abby irv_adam])
std_irv = std([irv_alice irv_abby irv_adam])

%% Table 1 Inspiratory Capacity

ic_alice = avg_vol_inspiration_alice + irv_alice
ic_abby = avg_vol_inspiration_abby + irv_abby
ic_adam = avg_vol_inspiration_adam + irv_adam

mean_ic = mean([ic_alice ic_abby ic_adam])
std_ic = std([ic_alice ic_abby ic_adam])

%% Table 1  ERV calculations
% alice ERV
[peaks, locs] = findpeaks(-alice_ex2_volume_smoothed); % find local minima
locs_alice = locs;
peaks_alice = peaks;
threshold = 0.45;  
minima_indices_alice = find(-alice_ex2_volume_smoothed(locs_alice) < threshold); 
locs_low_alice = locs_alice(minima_indices_alice);
figure
plot(time*dt, alice_ex2_volume_smoothed); % plot to verify correct peaks
hold on
plot(time(locs_low_alice)*dt, alice_ex2_volume_smoothed(locs_low_alice), "bo");
hold off

% do calculations
meanMediumMin_alice = mean(alice_ex2_volume_smoothed(locs_low_alice))
erv_alice = abs(min(alice_ex2_volume_smoothed) - meanMediumMin_alice)

% abby ERV
[peaks, locs] = findpeaks(-abby_ex2_volume_smoothed);
locs_abby = locs;
peaks_abby = peaks;
threshold = 0.5; % threshold changes since abby's data's minimum is much less
minima_indices_abby = find(-abby_ex2_volume_smoothed(locs_abby) < threshold); 
locs_low_abby = locs_abby(minima_indices_abby);
figure
plot(time*dt, abby_ex2_volume_smoothed);
hold on
plot(time(locs_low_abby)*dt, abby_ex2_volume_smoothed(locs_low_abby), "bo");
hold off
 
% do calculations
meanMediumMin_abby = mean(abby_ex2_volume_smoothed(locs_low_abby))
erv_abby = abs(min(abby_ex2_volume_smoothed) - meanMediumMin_abby)

% adam ERV
[peaks, locs] = findpeaks(-adam_ex2_volume_smoothed);
locs_adam = locs;
peaks_adam = peaks;
threshold = 0.45; % threshold doesn't seem to be applying; tried diff values
minima_indices_adam = find(-adam_ex2_volume_smoothed(locs_adam) < threshold); 
locs_low_adam = locs_adam(minima_indices_adam);
figure
plot(time*dt, adam_ex2_volume_smoothed);
hold on
plot(time(locs_low_adam)*dt, adam_ex2_volume_smoothed(locs_low_adam), "bo");
hold off

% do calculations
meanMediumMin_adam = mean(adam_ex2_volume_smoothed(locs_low_adam))
erv_adam = abs(min(adam_ex2_volume_smoothed) - meanMediumMin_adam)

% calculate mean + std for ERV across all group members
mean_erv = mean([erv_alice erv_abby erv_adam])
std_erv = std([erv_alice erv_abby erv_adam])


%% Table 1 EC calculations
ec_alice = avg_vol_inspiration_alice + erv_alice

ec_abby = avg_vol_inspiration_abby + erv_abby

ec_adam = avg_vol_inspiration_adam + erv_adam

mean_ec = mean([ec_alice ec_abby ec_adam])
std_ec = std([ec_alice ec_abby ec_adam])

%% Table 1 VC = IRV + ERV + Vt
vc_alice = irv_alice + erv_alice + avg_vol_inspiration_alice
vc_adam = irv_adam + erv_adam + avg_vol_inspiration_adam
vc_abby = irv_abby + erv_abby + avg_vol_inspiration_abby

mean_vc = mean([vc_alice vc_abby vc_adam])
std_vc = std([vc_alice vc_abby vc_adam])

%% Table 1 RV = predicted VC x 0.25 
rv_alice = vc_alice*0.25
rv_adam = vc_adam*0.25
rv_abby = vc_abby*0.25

mean_rv = mean([rv_alice rv_adam rv_abby])
std_rv = std([rv_alice rv_adam rv_abby])


%% Table 1 TLC =  VC + RV
tlc_alice = vc_alice + rv_alice
tlc_adam = vc_adam + rv_adam
tlc_abby = vc_abby + rv_abby

mean_tlc = mean([tlc_alice tlc_adam tlc_abby])
std_tlc = std([tlc_alice tlc_adam tlc_abby])

%% Table 1 FRC =  ERV + RV
frc_alice = erv_alice + rv_alice
frc_adam = erv_adam + rv_adam
frc_abby = erv_abby + rv_abby

mean_frc = mean([frc_alice frc_adam frc_abby])
std_frc = std([frc_alice frc_adam frc_abby])
