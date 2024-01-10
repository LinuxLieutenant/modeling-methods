% 
% EML3034C
% Project 9: 
% Due: 11-3-23

clc, clear, close all
format long

fprintf("=================================\n")
fprintf("\n")
fprintf("Project 9\n")
disp(datetime('today'))
fprintf("=================================\n")

%% definitions
%define constants (convert theta to radians)
tmax = 200;
u0 = [pi/6.25;0]; %recall that u contains theta and theta dot in a vector
printtime = 132.65;
L=4.5;
g=9.81;

%% Find the period ucing quad
%define k and period function T see project 7 copy paste.
theta = pi/6.25;
k = sin(0.5*(theta));
T = @(theta) ((4.*sqrt(L./g))./sqrt((1-(k).^2.*(sin(theta)).^2)));

%Define interval
a = 0;
b = pi/2;

%call quadgk
period = quadgk(T ,a, b, "RelTol", 1e-6);
fprintf("period = %.4e\n\n", period) %ensure period matches assignment

%% Using ODE45 to solve odes

tol = 1e-3;
options = odeset( "RelTol", tol, "AbsTol", tol); %set AbsTol and RelTol
[t3, thetas3] = ode45(@pendulum_drag, [0 tmax], u0, options);
%t3 holds the time steps across interval
%theta3 holds theta and thetadot values across interval

fprintf("using tol = %.1e\n", tol)

%use interp1 to evaluate theta dot at printtime
fprintf("Displacement is %.5f at t = %.3f\n", interp1(t3,thetas3(:,1),printtime), printtime)
fprintf("Velocity is     %.5f at t = %.3f\n", interp1(t3,thetas3(:,2),printtime), printtime)

%% Graph
%Define vector of points along pendalum period
T = 4.3237;
plot_T = t3;

%Displacement vs Time
figure(1)
plot(plot_T, thetas3(:,1)) %independent, dependent
hold on
plot(plot_T, thetas3(:,1), 'ro')
grid
title("Displacement from equilibrium vs t (tol 1e-3)")
legend ("theta(t)")
xlabel("Time (s)")
ylabel("Displacement")

%Displacement vs Time
figure(2)
plot(plot_T, thetas3(:,2)) %independent, dependent
hold on
plot(plot_T, thetas3(:,2), 'b')
grid
title("Mass Velocity vs t (tol 1e-3)")
legend ("Velocity")
xlabel("Time (s)")
ylabel("Theta dot(t)")

%Poincare Diagram
figure(3)
plot(thetas3(:,1), thetas3(:,2)) %independent, dependent
hold on
plot(thetas3(:,1), thetas3(:,2), 'b')
grid
title("Mass Velocity vs t (tol 1e-3)")
legend ("Trace")
xlabel("Theta (t)")
ylabel("Theta dot(t)")

%%Change tolerance and repeat
tol = 1e-6;
options = odeset("RelTol",tol,"AbsTol",tol); %set AbsTol and RelTol
[t6, thetas6] = ode45(@pendulum_drag, [0 tmax], u0, options);
%t3 holds the time steps across interval
%theta3 holds theta and thetadot values across interval

fprintf("using tol = %.1e\n", tol)

%use interp1 to evaluate theta dot at printtime
fprintf("Displacement is %.5f at t = %.3f\n", interp1(t6, thetas6(:,1), printtime), printtime)
fprintf("Velocity is     %.5f at t = %.3f\n", interp1(t6, thetas6(:,2), printtime), printtime)

plot_T = t6;
%Displacement vs Time
figure(4)
plot(plot_T, thetas6(:,1)) %independent, dependent
hold on
plot(plot_T, thetas6(:,1), 'ro')
grid
title("Displacement from equilibrium vs t (tol 1e-6)")
legend ("theta(t)")
xlabel("Time (s)")
ylabel("Displacement")

%Displacement vs Time
figure(5)
plot(plot_T, thetas6(:,2)) %independent, dependent
hold on
plot(plot_T, thetas6(:,2), 'b')
grid
title("Mass Velocity vs t (tol 1e-6)")
legend ("Velocity")
xlabel("Time (s)")
ylabel("Theta dot(t)")

%Poincare Diagram
figure(6)
plot(thetas6(:,1), thetas6(:,2)) %independent, dependent
hold on
plot(thetas6(:,1), thetas6(:,2), 'b')
grid
title("Mass Velocity vs t (tol 1e-6)")
legend ("Trace")
xlabel("Theta (t)")
ylabel("Theta dot(t)")

function [odes] = pendulum_drag(t, u)

%input u is a vector to couple the ODES and we never formally define it
%this routine takes t as an input to specify that theta and theta dot are 
%functions of time

%define constants
g = 9.81;
L = 4.5;
c = 0;
%c = 0.2; %comment out one value and run it then switch for the other value\
% of c for the quiz

%define our coupled ODEs in col vector
odes(1,1) = u(2);
odes(2,1) = -c*(u(2))*abs(u(2))-(g/L)*sin(u(1));

end