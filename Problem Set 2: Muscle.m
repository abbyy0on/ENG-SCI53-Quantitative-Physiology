%%
% File Name: pset2.m
% Author: Abby Yoon
% Created: Sep 8, 2023
% Description: Pset 2 Solutions

clear % Clears the workspace variables and Command Window
close all % Close all figure windows
clc % Clears your Command Window

%% 
  
% stimulate model neuron with a current pulse of 0.1 milisec duration and an
% amplitude of 100 uA/cm2; NOTE: I opened up the given MatLab files

pulse_height = 100;
pulse_width = 0.1;
T_final = 8;

[t,y] = run_hh_model(T_final,pulse_height,pulse_width);

%% plotting

GNa = 120;    % GNa = the max possible value of gNa (mS)
GK = 36;      % GK = the max possible value of gK

nNew  = 0.9494;
mNew  = 0.9953;
hNew = 0.0009;
n0  = y(:,4); % set up a new n0, using the initial conditions that were set for our state vectors in the other file
m0  = y(:,2);
h0 = y(:,3);
nTau = 1.2028;
mTau = 0.1577;
hTau = 1.0022;

n_t = nNew - (nNew - n0).*exp(-(t/nTau));
m_t = mNew - (mNew - m0).*exp(-(t/mTau));
h_t = hNew - (hNew - h0).*exp(-(t/hTau));


gNa = GNa.*m_t.^3.*h_t; % equations for conductances
gK = GK.*n_t.^4;

% membrane potential
subplot(3,1,1)
plot(t, y(:,1)) % this calls on all the rows in the first column
title("Membrane Potential of Stimulated Model Neuron")
ylabel("Membrane Potential (mV)")
xlabel("time(ms)")

% sodium and potassium conductances
subplot(3,1,2)
plot(t, gK, "k")
hold on 
plot(t, gNa, "r")
title("Potassium and Sodium Conductance")
ylabel("Conductance (ms/cm2)")
xlabel("time(ms)")
legend("Potassium", "Sodium")
hold off

% the m-gate, hgate, and n-gate open probabilities
subplot(3,1,3)
plot(t, m_t, "r")  
hold on 
plot(t, h_t, "k")
plot(t, n_t, "b")
title("M-Gate, H-Gate, and N-Gate Open Probabilities")
ylabel("Coefficient")
xlabel("time(ms)")
legend("M-Gate", "H-Gate", "N-Gate")
hold off

%% 1b 

% stimulate the hh_model 101 different times (use a "for" loop) with a
% current pulse of 0.1 msec duration & an amplitude varying from 0 to 100
% uA/cm2

figure, hold on
T_final = 10; % the figure on the p-set had a time scale of 10 ms
maxVoltage = [];
for c = 1:101
    pulse_height = c - 1;
    [t,y] = run_hh_model(T_final,pulse_height,pulse_width);
    plot(t, y(:,1))
    maxVoltage(end+1) = max(y(:,1)); % this adds a new value to an array
end
title("Membrane Voltage Traces")
ylabel("Membrane Potential (mV)")
xlabel("time(ms)")
hold off


%% 1c

% plot the relationship between the stimulating current pulse magnitude and
% peak/max membrane voltage for each voltage trace 

figure, hold on
pulse_height = (0:100);
plot(pulse_height, maxVoltage, ".-") 
title("All or None Behavior in Action Potential Firing")
ylabel("Peak Depolarization (mV)")
xlabel("Amplitude of a 0.1 ms Current Pulse (uA/cm^2)")

%% 1d
% choose a value where an action potential first happens - i.e., it seems
% to be 18 based on my observation while looking at the "maxVoltage" array

figure, hold on
[t,y] = run_hh_model(T_final,18,pulse_width);
plot(t, y(:,1), "r") 
[t,y] = run_hh_model(T_final,100,pulse_width);
plot(t, y(:,1), "k") 
title("Action Potentials from Different Current Pulses")
ylabel("Membrane Potential (mV)")
xlabel("Time (ms)")
legend("barely superthreshold", "fully superthreshold")

 
%% 3
% (3a) alter the code to change the frequency of excitation. converted to Hz
% later.
% ANSWER: The individual twitch can be seen at a frequency of 10 Hz. At 200
% Hz, you can reach tetany. Tetany is a smooth, straight line indicating a continuous contraction of the muscle. This is reached when dt=5

dt = 50;             % time period between twitches in ms
frequency = 1/dt;    % twitch frequency
numtwitches = 20;    % number of twitches to initiate
x = [0:0.1:40];  
twamp = zeros(1,2001);
twamps = zeros(1,2001);
 
% Generate time delayed twitches and sum
for i= 1:numtwitches;
y(i,:) = gampdf(x,3,1); %approximates a fast twitch response
start = round((i-1)*dt)+1;
twamp(i,start:(start+length(x)-1)) = y(i,:);
end
tetf = sum(twamp,1);
 
figure(5)
subplot(211)
plot(tetf,'r')
axis([0 500 0 1.1*max(max(tetf))])
xlabel('Time (ms)')
ylabel('Twitch Force (a.u.)')
legend('Fast Twitch','location','southeast')
 
% Now you generate time delayed slow twitches and sum them.
% Write your own code (or modify below) and plot your results in subplot(212)in blue
for j= 1:numtwitches;
y2(j,:) = gampdf(x,3,4); %approximates a slow twitch response
star = round((j-1)*dt)+1;
twamps(j,star:(star+length(x)-1)) = y2(j,:);
end
tets = sum(twamps,1);
  
subplot(212)
plot(tets,'b')
axis([0 500 0 1.1*max(max(tets))])
xlabel('Time (ms)')
ylabel('Twitch Force (a.u.)')
legend('Slow Twitch','location','southeast')

% (3b) change the number of twitches
% ANSWER: As you increase frequency,then you should reach tetany. If you decrease the number of twitches, then you can reach a point where you cannot see tetany, since tetany is the result of a sum of twitches. The force profile is reflected by this, as the graph shows no continuous contraction and instead individual twitches instead. 

% For other questions for 3, please see my PSET PDF Solutions

%% 5

% (5a) plot force-velocity relationship for the muscle
a = 20;
b = 0.2;
t0 = 100;
t = 0;
t = (0:100); % tension

% equation for velocity
v = (((t0 + a).*b) ./ (t+a))-b;
figure
plot(t, v);
xlabel('Force (tension)(N)')
ylabel('Velocity (m/s)')
title("Force-Velocity Curve")


% (5b) change a = 5 N to 30 N in steps of 5 N
t0 = 100;
maxV = 1;
figure
hold on 
for a = (5:5:30)
    b = (a*maxV)/t0;
    v = (((t0 + a).*b) ./ (t+a))-b;
    plot(t, v);
end
xlabel('Force (tension)(N)')
ylabel('Velocity (m/s)')
title("Force-Velocity Curve")
hold off



