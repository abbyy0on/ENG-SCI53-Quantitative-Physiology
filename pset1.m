%%
% File Name: pset.1
% Author: Abby Yoon
% Created: Sep 8, 2023
% Description: P-set 1 Answers

clear % Clears the workspace variables and Command Window
close all % Close all figure windows
clc % Clears your Command Window

%% 4a  

% Create equations for n(t), m(t), and h(t) and include all constants
% within a script

t = (0:8/200:8);    %  create a time vector  - the interval bounds are defined in 4c; this creates 201 data points
% assign it to variables
nNew  = 0.9494;
mNew  = 0.9953;
hNew = 0.0009;
nInitial = 0.3177;
mInitial = 0.0529;
hInitial = 0.5961;
nTau = 1.2028;
mTau = 0.1577;
hTau = 1.0022;

n = nNew - (nNew - nInitial)*exp(-(t/nTau));
m = mNew - (mNew - mInitial)*exp(-(t/mTau));
h = hNew - (hNew - hInitial)*exp(-(t/hTau));

%% 4b

% Create equations for gk(t) and gNa(t) 

gKMax = 36;
gNaMax = 120;
gK = (n.^4).* gKMax;    %  create an equation for gk(t)
gNa = (m.^3).* h.* gNaMax;  %  create an equation for gNa(t)

%% 4c
 
% Plot your functions for both gK(t) and gNa(t) as a function of t
plot(t, gK, "r")
title("Conductances after a voltage clamp from -65 to +23 mV")
ylabel("conductance (mS/cm2)")
xlabel("time(ms)")

hold on 
plot(t, gNa, "y")
hold off

%% 4d

% The graph of conductances differ from the red and yellow gK and gNa
% curves in that the K+ conductance curve does not reach a peak and decrease; instead mine seems to increase 
% gradually and reach a  bit of an asymptote within a similar time interval. This is likely
% because our conductances underwent a voltage clamp experiment where
% the membrane potential was infinitely clamped from -65 to +23 mV and
% cannot go beyond those boundaries. However, in our handout, the membrane
% potential was not specified to have any limits and were subjected to an action potential, in which
% the K+ channels were inactivated  and the conductance thus saw a gradual
% decrease because ions cannot pass through the membranes as easily. (We do not have voltage clamps in our body.) On
% the other hand, in the handout, a desired membrane potential was
% explicitly set for the experiment. Current was injected into the axon
% through an electrode to make the measured membrane potential equal to the
% desired potential. This current or flow of ions likely contributed to a
% higher K+ conductance that approaches 30 mV and adjusted what the
% conductance originally would have been.

%% 4e

% Create an expression for Vmem(t)

eK = -90;
eNa = 60;
vMem = ((gK .* eK) + (gNa .* eNa))./(gK + gNa);

figure
plot(t, vMem, "b")
xlabel("time(ms)")
ylabel("Vmem(t)")
title("Vmem(t)")

% NOTE: The Nernst potentials for each ion were obtained from pg. 8 of
% slide 22 in Sept 11 lecture. I forgot to do elementwise operations
% initially, and the graph looked like a flat line

% My graph of Vmem(t) differs from what I'd expect for an action potential,
% because an action potential occurs when  the cells are excited by a
% stimulus that induces depolarization, hyperpolarization, etc. This
% influences the resting membrane potential and changes it into different
% values as time passes. On the other hand, Vmem(t) is the membrane
% potential when cells are not excited and are at rest. So it makes sense
% that they are different.

%% 5

% Declare a time vector from -30 to 30, with increments of 0.001
t = (-30:0.001:30);

% Plot a sinc function using the time vector, in black 
figure
y = sinc(t);
plot(t, y, "k");
title("Sinc Function");
ylabel("sinc(t)");
xlabel("Time (t)");
legend("sinc(t)");

% Plot the derivative of sinc function using the time vector, in black;
% "diff" is one code to use, but "gradient()" works as well
figure
diffY = gradient(y);
plot(t, diffY, "k");
title("Derivative of Sinc Function");
ylabel("d(sic(t))/dt");
xlabel("Time (t)");
legend("d(sinc(t))/dt");

%  Third figure:
% Plot the positive portions of this derivative with green points
% Plot the negative portions of this derivative with red points
% Mark the zero-crossing points of the derivative with yellow squares
signDY  =  sign(diffY);   

figure, hold on     % NOTE: You don't need to do 'hold on' so many times! It's like an on-off switch. Also, make sure to do "figure" because it creates a blank canvas
plot(t(signDY == 1), diffY(signDY == 1), "g.");
plot(t(signDY == 0), diffY(signDY == 0), "ys");
plot(t(signDY == -1), diffY(signDY == -1), "r.");
title("Derivative of Sinc Function");
ylabel("d(sinc(t))/dt");
xlabel("Time (t)");
legend("Positive derivative", "Zero-crossing point derivative", "Negative derivative");
hold off

% Plot the original sinc function. In one figure:
% Plot the portions of the original function where its derivative is positive with green points
% Plot the portions of the original function where its derivative is negative with red points
% Mark the peaks and troughs of the original function with yellow squares
 

figure, hold on    
plot(t(signDY == 1), y(signDY == 1), "g.");
plot(t(signDY == 0), y(signDY == 0), "ys");
plot(t(signDY == -1), y(signDY == -1), "r.");
title("Sinc Function");
ylabel("sinc(t)");
xlabel("Time (t)");
legend("Positive derivative", "Peaks/troughs", "Negative derivative");
hold off