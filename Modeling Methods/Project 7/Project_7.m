% 
% EML3034C
% Project 7: 
% Due: 11-3-23

clc, clear, close all
format long

fprintf("=================================\n")
fprintf("\n")
fprintf("Project 7\n")
disp(datetime('today'))
fprintf("=================================\n")

%definitions
%define constants (convert theta to radians)
theta=28.8*(pi/180);
k=sin(theta./2);
L=14.5;
g=9.81;

%define pendulum function and test function

T = @(theta) ((4.*sqrt(L./g))./sqrt((1-(k).^2.*(sin(theta)).^2)));
test = @(x) ((12.5)*x.^9+(3.2).*x.^4);

%Test Integral

fprintf("Test function quadrature\n");
fprintf("------------------------\n");

%define a, b
a=1.2;
b=2.5;

%call GaussQuad
val = GaussQuad(test,a,b,5);
%print values for test run
fprintf("a = %.3e \n", a)
fprintf("b = %.3e \n", b)
fprintf("ans = %.10e \n", val)

%use MATLAB function "quad" to solve and compare
%print MATLAB answers
fprintf("MATLAB = %.10e \n", quad(test, a, b))

%find exact test integral
%exact antiderivative:
antiderivative = @(x) (((12.5)*x^10)/10)+(((3.2)*x^5)/5);
%find exact integral across interval:
integral = ((((12.5)*b^10)/10)+(((3.2)*b^5)/5)) - ((((12.5)*a^10)/10)+(((3.2)*a^5)/5));
%print out exact test results
fprintf("Exact answer = %.10e \n", integral) %from antiderivative
fprintf("Absolute error = %.10e \n", abs(integral-val)) %error exact answer and GaussQuad
fprintf("Relative (percent) error = %.10e \n", (abs(integral-val)/integral)) %abs error / exact answer

%Gauss-Legendre 5-point Quadrature in 1 interval
fprintf("One interval: 0 to pi/2 \n")
fprintf("------------------------\n")

%define a and b 
a = 0;
b = pi/2;

%call GaussQuad
val = GaussQuad(T,a,b,5);

fprintf("a = %.3e \n", a)
fprintf("b = %.3e \n", b)
fprintf("ans = %.10e \n\n", val)


%Gauss-Legendre 5-point Quadrature in 2 intervals
%first half

fprintf("Two intervals: 0 to pi/4 and pi/4 to pi/2\n")
fprintf("-----------------------------------------\n")

%define a and b 
a = 0;
b = pi/4;
%call GaussQuad
val1 = GaussQuad(T,a,b,5);
%solve using MATLAB
val1M = quad(T,a,b);
fprintf("a = %.3e \n", a)
fprintf("b = %.3e \n", b)
fprintf("ans = %.10e \n\n", val1)

%second half

%define a and b 
a = pi/4;
b = pi/2;
%call GaussQuad
val2 = GaussQuad(T,a,b,5);
%solve using MATLAB
val2M = quad(T,a,b);
fprintf("a = %.3e \n", a)
fprintf("b = %.3e \n", b)
fprintf("ans = %.10e \n\n", val2)

%sum the two halves to find overall answer across interval:
fprintf("Our answer = %.10e seconds.\n", val1+val2)
fprintf("MATLAB's answer = %.10e seconds. \n", val1M+val2M)

function [val] = GaussQuad(f, a, b, num_points)
  
  if num_points==2
    % 2 point Gauss
    fprintf('using 2 points Gauss Quad\n')
    % weights
    w = ones(2,1); % same as w = [1; 1]
    % gauss points
    xi = zeros(2,1);
    xi(1) = 1/sqrt(3);
    xi(2) = -xi(1);
  elseif num_points==5
    % 5 point Gauss
    fprintf('using 5 points Gauss Quad\n')
    % weights (Based on the provided table in the question)
    w = zeros(5,1);
    w(1) = 0.236926885056189;
    w(2) = 0.478628670499366;
    w(3) = 0.568888888888889;
    w(4) = 0.478628670499366;
    w(5) = 0.236926885056189;
    % gauss points
    xi = zeros(5,1);
    xi(1) = 0.906179845938664;
    xi(2) = 0.538469310105683;
    xi(3) = 0.000000000000000;
    xi(4) = -0.538469310105683;
    xi(5) = -0.906179845938664;
  else
    fprintf('did not enter 2 or 5\n')
  end
    %calculate jacobian
    J=(b-a)/2;
    %map xi to domain of integration
    x=J*xi+(b+a)/2;
    %perform quadrature
    val=J*sum(w .* f(x));
    
  end