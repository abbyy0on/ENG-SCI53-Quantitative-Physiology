%% ENG-SCI 53 Lab 3
% Filename: es53_lab4_laf.m
% Author: Abby Yoon, Lydia Shimelis, Fredrik Willumsen Haug
% Created: October 14th
% Description: Blood Pressure and Blood Flow

%%
%% Figure 1
% load data and define samplerate
fredrik = load('fredrik_ex2.mat');
sr = fredrik.samplerate(1);
% fredrik blood pressure and pulse data
fredrik_blood_pressure = fredrik.data(fredrik.datastart(1, 1) :
fredrik.dataend(1, 1));
fredrik_pulse = fredrik.data(fredrik.datastart(2, 1) : fredrik.dataend(2,
1));
timef = (0 : length(fredrik_blood_pressure)-1);
% plots within a 1x2 subplot
figure(1);
hold on
subplot(2, 1, 1);
plot(timef / sr, fredrik_blood_pressure);
xlabel("Time (s)");
ylabel("Blood Pressure (mmHG)");
hold off
figure(1);
hold on
subplot(2, 1, 2);
plot(timef / sr, fredrik_pulse);
xlabel("Time (s)");
ylabel("Pulse Wave Magnitude");
hold off

%%
%% Figure 2
% load data and define samplerate
abby = load('abby_ex2.mat');
sr = abby.samplerate(1);
% specific blood pressure and pulse data based on channels and blocks
abby_blood_pressure = abby.data(abby.datastart(1, 1) : abby.dataend(1,
1));
abby_pulse = abby.data(abby.datastart(2, 1) : abby.dataend(2, 1));
timea = (0 : length(abby_blood_pressure)-1);
% plots within a 1x2 subplot
figure(2);
hold on
subplot(2, 1, 1);
plot(timea / sr, abby_blood_pressure);
xlabel("Time (s)");
ylabel("Blood Pressure (mmHG)");
hold off
figure(2);
hold on
subplot(2, 1, 2);
plot(timea / sr, abby_pulse);
xlabel("Time (s)");
ylabel("Pulse Wave Magnitude");
hold off

%%
%% Figure 3
% load data and define samplerate
lydia = load('lydia_ex2.mat');
sr = lydia.samplerate(1);
% specific blood pressure and pulse data based on channels and blocks
lydia_blood_pressure = lydia.data(lydia.datastart(1, 1) : lydia.dataend(1,
1));
lydia_pulse = lydia.data(lydia.datastart(2, 1) : lydia.dataend(2, 1));
timel = (0 : length(lydia_blood_pressure)-1);
% plots within a 1x2 subplot
figure(3);
hold on
subplot(2, 1, 1);
plot(timel / sr, lydia_blood_pressure );
xlabel("Time (s)");
ylabel("Blood Pressure (mmHG)");
hold off
figure(3);
hold on
subplot(2, 1, 2);
plot(timel / sr, lydia_pulse);
xlabel("Time (s)");
ylabel("Pulse Wave Magnitude");
hold off

%%
%% Blood Pressure Extraction Based on Korotkoff sound start / end
% section load in case it's helpful
% fredrik = load('fredrik_ex2.mat');
% abby = load('abby_ex2.mat');
% lydia = load('lydia_ex2.mat');
% fredrik blood pressure and mean arterial pressure
fredrik_systolic = fredrik.data(fredrik.com(1, 3));
fredrik_diastolic = fredrik.data(fredrik.com(2, 3));
fredrik_bp_value = fredrik_systolic / fredrik_diastolic;
fredrik_map = ((2/3) * fredrik_diastolic) + ((1/3) * fredrik_systolic);
% abby blood pressure and mean arterial pressure
abby_systolic = abby.data(abby.com(1, 3));
abby_diastolic = abby.data(abby.com(2, 3));
abby_bp_value = abby_systolic / abby_diastolic;
abby_map = ((2/3) * abby_diastolic) + ((1/3) * abby_systolic);
% lydia blood pressure and mean arterial pressure
lydia_systolic = abby.data(lydia.com(1, 3));
lydia_diastolic = fredrik.data(lydia.com(2, 3));
lydia_bp_value = lydia_systolic / lydia_diastolic;
lydia_map = ((2/3) * lydia_diastolic) + ((1/3) * lydia_systolic);
% mean and std of digital group mean arterial pressures
group_map = [fredrik_map, abby_map, lydia_map];
mean_group_map = mean(group_map);
std_group_map = std(group_map);

%%
%% Figure 4
% section load in case it's helpful
% fredrik = load('fredrik_ex2.mat');
% abby = load('abby_ex2.mat');
% lydia = load('lydia_ex2.mat');
% samplerate declaration
S = fredrik.samplerate(1);
% NB! Plot in report includes all data because it resembles the
% type of plot in the lab handout swithch the comments for the two
% blocks below to get either full data plot or only peak relevant for
% HR and discard HR calculations if using first block for general plot
% fredrik_pulse = fredrik.data(fredrik.datastart(2, 1): fredrik.dataend(2,
1));
% abby_pulse = abby.data(abby.datastart(2, 1) : abby.dataend(2, 1));
% lydia_pulse = lydia.data(lydia.datastart(2, 1) : lydia.dataend(2, 1));
% pulse for each subject extracted when peaks are apparent due to present
% pulse
fredrik_pulse = fredrik.data(fredrik.datastart(2, 1) + 7*S:
fredrik.datastart(2, 1) + 19*S);
abby_pulse = abby.data(abby.datastart(2, 1) : abby.datastart(2, 1) +
13*S);
lydia_pulse = lydia.data(lydia.datastart(2, 1) + 32*S : lydia.datastart(2,
1)+45*S);
% fredrik peaks
e20_1 = fredrik_pulse / max(fredrik_pulse);
de20_1 = diff(e20_1); % derivative of this vector
threshold = 0; % peak treshold
ind_1 =
find((e20_1(2:end-1)>threshold)&(de20_1(1:end-1)>0)&(de20_1(2:end)<0))+1;
% the indices of all peaks in the first 20 second segment
% abby peaks
e20_abby = abby_pulse / max(abby_pulse);
de20_abby = diff(e20_abby); % derivative of this vector
threshold = 0; % peak treshold
ind_1 =
find((e20_abby(2:end-1)>threshold)&(de20_abby(1:end-1)>0)&(de20_abby(2:end
)<0))+1; % the indices of all peaks in the first 20 second segment
% abby peaks
e20_lydia = lydia_pulse / max(lydia_pulse);
de20_lydia = diff(e20_lydia); % derivative of this vector
threshold = 0; % peak treshold
ind_1 =
find((e20_lydia(2:end-1)>threshold)&(de20_lydia(1:end-1)>0)&(de20_lydia(2:
end)<0))+1; % the indices of all peaks in the first 20 second segment
% fredrik plot
time_1 = (1:length(e20_1))/S;
figure(4);
subplot(3, 1, 1);
plot(time_1, e20_1);
hold on
xlabel('time (s)');
ylabel('PWM Fredrik');
axis([0 10 -1.1 1.1]);
MPH = threshold;
MPD = S/2;
[~,locs_Rwave_s1] =
findpeaks(e20_1,'MinPeakHeight',MPH,'MinPeakDistance',MPD);
plot(time_1(locs_Rwave_s1), e20_1(locs_Rwave_s1),'r*');
legend('Pulse Wave', 'Peak');
hold off
% abby plot
timea = (1:length(e20_abby))/S;
figure(4);
subplot(3, 1, 2);
plot(timea, e20_abby);
hold on
xlabel('time (s)');
ylabel('PWM Abby');
axis([0 10 -1.1 1.1]);
MPH = threshold;
MPD = S/2;
[~,locs_Rwave_s2] =
findpeaks(e20_abby,'MinPeakHeight',MPH,'MinPeakDistance',MPD);
plot(timea(locs_Rwave_s2), e20_abby(locs_Rwave_s2),'r*');
legend('Pulse Wave', 'Peak');
hold off
% lydia plot
timel = (1:length(e20_lydia))/S;
figure(4);
subplot(3, 1, 3);
plot(timel, e20_lydia);
hold on
xlabel('Time (s)');
ylabel('PWM Lydia');
axis([0 10 -1.1 1.1]);
MPH = threshold;
MPD = S/2;
[~,locs_Rwave_s3] =
findpeaks(e20_lydia,'MinPeakHeight',MPH,'MinPeakDistance',MPD);
plot(timel(locs_Rwave_s3), e20_lydia(locs_Rwave_s3),'r*');
legend('Pulse Wave', 'Peak');
hold off
% heart rate calculations
% fredrik
peak_interval_f = diff(locs_Rwave_s1) / S;
mean_peak_interval_f = mean(peak_interval_f);
mean_heart_rate_fredrik = 60 / mean_peak_interval_f
% abby
peak_interval_a = diff(locs_Rwave_s2) / S;
mean_peak_interval_a = mean(peak_interval_a);
mean_heart_rate_abby = 60 / mean_peak_interval_a
% lydia
peak_interval_l = diff(locs_Rwave_s3) / S;
mean_peak_interval_l = mean(peak_interval_l);
mean_heart_rate_lydia = 60 / mean_peak_interval_l

%%
%% Figure 7
expected_mean = 8.5; % value from table L or R of mean flow for 5-20 yo
subjects_data = [10.1610887, 7.4220126, 6.7858401];
% calculate the standard deviations
expected_std = 1.3;
subject_std = std(subjects_data);
% make the figure
figure(7)
y = [expected_mean, mean(subjects_data)];
x = ["Expected", "Experimental Average"];
% make a bar plot
bar1 = bar(y);
% add error bars using the standard deviations
hold on;
errorbar(1:2, y, [expected_std, subject_std], 'k', 'linestyle', 'none');
xticklabels(x);
% label plot
ylabel('Flow (mL/min)')
xlabel('Subject')
hold off;

%%
%% Table 4: Calculating Cardiac Output and Stroke Volume
% heart rates as calculated from before section
hr_f = 87;
hr_a = 75;
hr_l = 84;
% calculate CCA flow by mutiplying velocity (mL/s) and cross-sectional
area
% of the carotid atery, calculated in previous section in report
%cca_f = 48.2548; % v_f * a_f
%cca_a = 26.4581; %v_a * a_a
%cca_l = 17.671458; % v_l * a_l
cca_f = 10.16 * 60; % v_f * a_f in mL/min
cca_a = 7.42 * 60; %v_a * a_a in mL/min
cca_l = 6.79 * 60 % v_l * a_l in mL/min
% 10 percent of flow to the brain through carotid arteries (in OH, Dr.
% Moyer said that it would be 5% in the artery because it splits so
% multiply by 20)
co_f = (20 * cca_f)
co_a = (20 * cca_a)
co_l = (20 * cca_l)
co_mean = mean([co_f co_a co_l]) % calculate mean and std
co_std = std([co_f co_a co_l])
% calculate stroke volume
sv_f = co_f/hr_f
sv_a = co_a/hr_a
sv_l = co_l/hr_l
sv_mean = mean([sv_f, sv_a, sv_l])
sv_std = std([sv_f, sv_a, sv_l])
% calculate EDV using stroke volume and estimate that EF = .6 (60%)
edv_f = sv_f/.6
edv_a = sv_a/.6
edv_l = sv_l/.6
% calculate esv using EDV and SV
esv_f = edv_f - sv_f
esv_a = edv_a - sv_a
esv_l = edv_l - sv_l

%%
%% Figure 8
figure(8)
% fredrik data
pv_f = [93.43 233.56 233.56 93.43; 0 0 75 126]; % used bp from analog and
ESV and EDV from above
volume_f = [pv_f(1,:) pv_f(1,1)]; % close loop with first point
pressure_f = [pv_f(2,:) pv_f(2,1)]; % close loop with first point
% abby data
pv_a = [79.15 197.87 197.87 79.15; 0 0 67 90];
volume_a = [pv_a(1,:) pv_a(1,1)];
pressure_a = [pv_a(2,:) pv_a(2,1)];
% lydia data
pv_l = [64.67 161.67 161.67 64.67; 0 0 67 104];
volume_l = [pv_l(1,:) pv_l(1,1)];
pressure_l = [pv_l(2,:) pv_l(2,1)];
% plot!!
plot(volume_f, pressure_f,'LineWidth', 2);
hold on;
plot(volume_a, pressure_a,'LineWidth', 2);
plot(volume_l, pressure_l,'LineWidth', 2);
xlabel("Volume (mL)");
ylabel("Pressure (mmHg)");
% legend
legend("Fredrik PV Loop", "Abby PV Loop", "Lydia PV Loop");
hold off;

%%
%% Figure 9: Calculating Stroke Work
% input variables for diasolic and systoic blood pressures
sys_f = 126; % used analog since that seemed more reasonable
dia_f = 75;
sys_a = 90;
dia_a = 67;
sys_l = 104;
dia_l = 67;
sys_e = 120;
dia_e = 80;
% calculate MAP given diastolic and systolic bp (Mean Arterial Presure =
% 1/3(systolic bp) + 2/3(diastolic bp))
map_e = (sys_e/3)+(2*dia_e/3)
map_f = (sys_f/3)+(2*dia_f/3)
map_a = (sys_a/3)+(2*dia_a/3)
map_l = (sys_l/3)+(2*dia_l/3)
% calculate stroke work for each member using stroke volume and MAP
(Stroke
% work = MAP (mmHg) * SV (mL))
sw_e = (80 * map_e) % calculated from class where 80 is the stroke volume
sw_f = sv_f * map_f
sw_a = sv_a * map_a
sw_l = sv_l * map_l
% convert to joules (stroke work (mmHg*mL) * 133 (N/m^2) * 10^(-6)m^3)/mL)
sw_e_j = sw_e * .133 * 10 ^(-3)
sw_f_j = sw_f * .133 * 10 ^(-3)
sw_a_j = sw_a * .133 * 10 ^(-3)
sw_l_j = sw_l * .133 * 10 ^(-3)
y = [sw_e_j sw_f_j sw_a_j sw_l_j];
figure (9)
bar1 = bar(y);
% Set the x-axis labels
x = ["Expected" "Fredrik" "Abby" "Lydia"];
xticks(1:numel(x)); % set the x-ticks at positions 1 and 2
xticklabels(x); % set the x-tick labels to the categories in 'x'
ylabel("Stroke Work (J)")
xlabel("Subject")
% calculate total work in one minute (formula: Heart rate (bpm) * SW J)
total_work_e = sw_e_j * 80 % 80 is middle of 60-100 bpm although closer to
60 at rest for "normal adult"
total_work_f = sw_f_j * hr_f
total_work_a = sw_a_j * hr_a
total_work_l = sw_l_j * hr_l
%%