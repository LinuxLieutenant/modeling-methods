% Project 8 - ODE Solver for Euler
% EML 3034C
% Due 11-10-23
% 
clc, clear, close all

fprintf("EML3034C Project 8\n\n\n");

%% Parameters
% Define constants
T0 = 298;
p = 16600;
cp = 160;
h = 125;
D = 1.5*10^-3;
B = 0.35;
eps = 0.925;
sig = 5.68 * 10^-8;

A = 4*pi*(D/2)^2;
V = (4/3)*pi*(D/2)^3;

%% Define time constant 
tau = 1/(h*A/(p*cp*V));
% Choosing timestep for convergence:
% for stability dt <= tau
% For strict undershoot, dt <= tau/2
fprintf("Time Constant: dt <= %.3f\n", tau)
fprintf("For Strict Undershoot: dt <= %.3f\n", tau/2)

%% Define ODE
% Define non-constant parameters
T_inf = @(t) T0 + B*t; % Function of time
hr = @(t, T) eps*sig*(T^2 + (T_inf(t))^2)*(T + T_inf(t)); % Function of two variables
T_rad= @ (t,T) ((A*(h+hr(t,T)))/(p*cp*V));
dTdt = @ (t,T) T_rad(t,T)*(T_inf(t))-(T_rad(t,T)*T);
%f = @(T, t) (((A*(h+hr(t,T)))/(p*cp*V)) * T_inf(t)) - (((A*(h+hr(t,T)))/(p*cp*V)) * T);
%dTdt = f; % Function of two variables


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% using Euler's Method

dt = [1, 5, 10];
tf = 200;
print_time = [50, 150, 200]; % Fill vector with desired time values
tol = [10^-3, 10^-6]; % Fill vector with tolerances for ode45

for i =1:length(dt) % Loop through time steps

    fprintf('using dt = %.4f\n',dt(i));
    [time, sol] = euler_solver(dTdt, dt(i), tf, T0);

    % interpolate solution to get value at time point we want
    for j = 1:length(print_time) % Loop through print_time

        interp_sol = interp1(time, sol, print_time(j));
        fprintf('Euler: %.4f at t = %.4f\n', interp_sol, print_time(j));
    end

    % plot the solutions
    figure
    plot(time, sol, 'b', 'linewidth', 3)
    title(sprintf("Temperature vs. Time for dt = %.3f\n",dt(i)));
    fprintf('\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% using adaptive Runge-Kutta 4th order method

for i = 1:length(tol) % Loop through tol vector
    
    fprintf('using tol = %.4e\n',tol(i));
    opt = odeset("RelTol", tol(i), "AbsTol", tol(i));
    [time_MATLAB, sol_MATLAB] = ode45(dTdt, [0, tf], T0, opt);
    
    % interpolate solution to get value at time point we want
    for j = 1:length(print_time) % Loop through print time
        interp_sol = interp1(time_MATLAB, sol_MATLAB, print_time(j));
        fprintf("ODE45: %.4f at t = %.4f\n", interp_sol, print_time(j));
    end
    
    % plot the solutions
    figure
    plot(time_MATLAB,sol_MATLAB,'b','linewidth',3)
    title(sprintf("Temperature vs. Time for tol = %.3e\n",tol(i)));
    fprintf('\n');
end

function [t, y] = euler_solver(f, dt, tf, y0)

    % Find the number of steps, N
    % Round up in MATLAB using ceil() always rounds up to nearest integer
    N = ceil(tf/dt);
    % Discretize time vector, t, using number of steps, N
    t=0:dt:dt*N;
    % Define y vector using initial condition
    y(1)=y0;
    % Loop through time and perform Euler's method
    for i=1:N
        y(i+1)=y(i)+dt*f(t(i),y(i));
    end
end