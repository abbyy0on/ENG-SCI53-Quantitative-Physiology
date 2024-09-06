%%
% File Name: pset4.m
% Author: Abby Yoon
% Created: October 15, 2023
% Description: Pset 4 Solution

clear % Clears the workspace variables and Command Window
close all % Close all figure windows
clc % Clears your Command Window

load('PVloop.mat');

%% 4a
volume = PVloop(1,:); % extract first row
pressure = PVloop(2,:); % extract second row

% make a complete loop
volume(end+1) = 60; % 60 mL is the volume that the loop should close at 
pressure(end+1) = 3.5; % 3.5 mmHg is the pressure that the loop should close at

plot(volume,pressure, "k");   
xlim([0 160]); % adjust x axis scale
ylim([0 160]); % adjust y axis scale
title('Pressure-Volume Loop');
xlabel('Volume (mL)');
ylabel('Pressure (mmHg)');

%% 4b
% NOTE: stroke work  = stroke volume * MAP

% SV = EDV - ESV
edv = max(volume);
esv = min(volume);
sv = edv - esv;

% MAP = 1/3(systolic pressure) + 2/3(diastolic pressure)
systolic = max(pressure);   % find systolic

v_index = find(volume == edv);  % find diastolic
pressure_at_v_index = pressure(v_index);
diastolic = pressure_at_v_index(2); 

map = (1/3)*systolic + (2/3)*diastolic;     % find MAP

% calculate stroke work
strokeWork = map*sv;
% convert into Joules (1 mmHg·mL = 1.3332*10^-4 J)
strokework_in_joules = strokeWork*(1.3332*10^-4)

%% 4c
% area inside PV loop
area = polyarea(volume, pressure);
% convert into Joules (1 mmHg·mL = 1.3332*10^-4 J)
area_in_joules = area*(1.3332*10^-4)


