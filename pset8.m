%%
% File Name: pset8.m
% Author: Abby Yoon

clear % Clears the workspace variables and Command Window
close all % Close all figure windows
clc % Clears your Command Window

%% a
% Define the range of x values
t_a = 0:.1:10;
% Compute e^(-x) for each x value
function_a = exp(-t_a);
asterisk_a = find(t_a == 1);
% Plot
figure
plot(t_a, function_a);
hold on
plot(t_a(asterisk_a), function_a(asterisk_a), 'r*');
hold off
grid on
xlabel("t")
ylabel('$e^{-t}$', 'Interpreter', 'latex')
legend('$e^{-t}$', 't=1', 'Interpreter', 'latex')

%% b
t_b = 0:.1:10;
function_b = (1) - (exp(-t_b));
asterisk_b = find(t_b == 1);
figure;
plot(t_b, function_b);
hold on
plot(t_b(asterisk_b), function_b(asterisk_b), 'r*');
hold off
grid on
xlabel("t")
ylabel('$1 - e^{-t}$', 'Interpreter', 'latex')
legend('$1 - e^{-t}$', 't=1', 'Interpreter', 'latex')

 
%% c
t_c = 0:.1:20;
a_0 = 2;
tau = 3;
A = a_0*exp(-t_c/tau);
figure;
plot(t_c, A);
hold on
hold off
grid on
xlabel("t")
ylabel('A')

 
%% d  

% log(X) returns ln(X) of each element
% time at A = 1/2 (A_0)
t_1 = log(0.5) * -tau
% time at A = .37 (A_0)
t_1 = log(0.37) * -tau

%% e
t_c = 0:.1:40;
a_0 = 2;
tau = 15;
A_2 = a_0*exp(-t_c/tau);
figure;
plot(t_c, A_2);
hold on
hold off
grid on
xlabel("t")
ylabel('A')

