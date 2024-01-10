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