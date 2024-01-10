% Project 4
% 
% EML3034C
% Due 10/6/2023
clc, clear, close all;
format long;

fprintf("=================\n")
fprintf("\n")
fprintf("Project 4: Gauss Seidel SOR\n")
disp(date())
fprintf("==============\n")

%A = 4*4 matrix ;
%b = 4*1 vector;
A= dlmread("fin_data_A_matrix_FA_2023.txt"); %make sure file name matches
b= dlmread("fin_data_b_vector_FA_2023.txt"); %File must be in working directory

% initial guess
guess = zeros(length(b),1);
guess(:) = 175;

% residual and iterative tolerance
tol = 1e-4;
max_iters = 1e+5;

% SOR parameter
omega = 1;

% run Gauss-Seidel with SOR
[x,res,conv,num_iters] = gs_sor(A, b, guess, tol, omega, max_iters);

fprintf("using omega = %.2f", omega)
fprintf('\nnum iters to converge = %i\n',num_iters)

%for testing only!
%fprintf('x = \n')
%disp(x)

%print for quiz
fprintf('T(20) is %f \n', x(20))
fprintf('T(58) is %f \n', x(58))
fprintf('T(80) is %f \n', x(80))
fprintf('T(123) is %f \n', x(123))

fprintf('\nres = %.6e\n', res(end))
fprintf('conv = %.6e\n\n', conv(end))

%compare test vector to MATLAB
x_MATLAB = linsolve(A,b);
%fprintf('x_MATLAB = \n')
%disp(x_MATLAB)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("===================================\n")

% run GS SOR for different omegas and save number of iterations
%  to convergence for each omega
omegas = 1:0.1:1.9;  % remember, start1.1 : step : stop1.9

for i = 1:length(omegas)

       %run GS_SOR
       [x, res, conv, n] = gs_sor(A, b, guess, tol, omegas(i), max_iters); %use a diff omega each time

       %make a vector to store the number of iterations per omega
       sor_iters(i) = n;
       %print out current omega and num iters to solve
       fprintf("Omega: %.1f \t Iterations: %i\n", omegas(i), sor_iters(i))
end
%fprintf('x=%.2f\n', x(123))

% plot omega vs number of iterations to convergence
figure(1)  % number your figures as below so you can always reference it later
plot(omegas, sor_iters) %x,y
%make a title and labels

title("Finding the Optimal Omega")
xlabel("Relaxation Parameters -> Omega")
ylabel("Numbers of Iterations to Converge")

% plot residual norm at each iteration for different omegas
% we could automate this but since only a few lines.. 
%  probably faster and simpler to copy/paste
%  different story if we wanted to do this many times..
 
% use "semilogy" instead of "plot" to plot y axis using log10 scale
% feel free to use "plot" instead of "semilogy" so your plots look like
%  in lecture (decaying exponential)
figure(2)
hold on
set(gca, 'Xscale', 'log', 'Yscale', 'log')

for i = 1.1:0.1:1.6 %range of omegas 1.1 to 1.6
    %call GS SOR
    [x, res, conv, n] = gs_sor(A, b, guess, tol, i, max_iters);
    semilogy(1:n,res)
    %plot = (%iters, blank)
end

grid on
%include title and axis labels
title("Residuals at each iteration for various omegas")
xlabel("Iteration")
ylabel("Residual Infinity Norm")
legend("1.1", "1.2", "1.3", "1.4", "1.5", "1.6")

function [x,res,conv,n] = gs_sor(A,b,guess, tol, omega, max_iters)

    x = guess;
    n = 1; %Iteration counter
    while n < max_iters
        x_old = x; %store old solutions

        for i = 1:length(b)
            %INSERT GS formula
            %x(i) = (1 / A(i, i)) * (b(i) - A(i, 1:i-1) * x(1:i-1) - A(i, i+1:end) * x_old(i+1:end));
            x(i)= (b(i)-A(i,:)*x(:)+A(i,i)*x(i))/(A(i,i));
            %Jacobi
            % x(i) = (1 / A(i, i)) * (b(i) - A(i, :) * x_old + A(i, i) * x_old(i));
            %INSERT SOR Formula
            x(i) = omega * (x(i)) + (1-omega) * x_old(i);
            %x(i)=x_old(i)+(omega*(x(i)-x_old(i)));

        end

    %calculate residual and store inf norm
    %Only want to store one res value (norm) per iteration
    res(n) = norm(A*x-b, inf);
    %Find conv norm
    conv(n) = norm(abs(x-x_old),inf);

    %test fpr convergence
    if res(n) < tol && conv(n) < tol
        break
    end

    %check for divergance
    if res(n) > 1e15
        fprintf("\nSolution has diverged!!!\n")
        break
    end

    n = n+1;
    
    end

end
