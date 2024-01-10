%
%EML 3034C Modeling Methods
%Project 5
%Due 10-13-2023

clc, clear, close all
format long

fprintf("======================================\n")
fprintf("Project 5: NLNR\n")
fprintf("\n")
disp(datetime("today"))
fprintf("======================================\n")

%Initial Guess
x = [1;1;1;1];
%x = [ ; ; ; ]
%x = [ , , , ]

N = length(x);
maxiters = 100;
tol = 10^-5;

res = zeros(N,1); %allocate space
resnorm = zeros(N,1);
J = zeros(N);
b = zeros(N,1);

for k = 1:maxiters %iteration counter

    fprintf("Iteration number %i\n", k) %fill with iteration number

    %allocate space
    %J = zeros(N);
    %b = zeros(N,1);

    for i = 1:N %loop through rows of Jacobian J and rhs vector b
        for j = 1:N %loop through columns of J
            J(i, j) =  Jacobian(i, j, x);%fill Jacobian here
        end

        b(i,1) = springfuncs(i, x); %fill rhs with F function values

    end

    x_old = x;
    dx = J\(-b); %only one negative for b
    %dx = linsolve(J, -b);
    %allows us to solve for dx (change in x aka conv)
    x = x_old + dx; %update answer

    %res = zeros(N,1); %allocate space
    %resnorm = zeros(N,1);

    for i = 1:N 
        res(i) = springfuncs(i, x); %residual = f(x) just like previous NR
    end

    %store res norm for every iteration
    resnorm(k) = norm(res, inf);

    %print results for current iteration
    fprintf("Residual Norm %.8e\n\n", resnorm(k))

    fprintf("Jacobian Matrix at current iteration\n")
    fprintf("---------------------------------------\n")
    format long g
    disp(J); %display Jacobian Matrix

    fprintf("----------------------------------------\n")
    fprintf("Approximated Vector Solution\n")
    fprintf("%.8f\n", x(:,1)) %print solution vector
    fprintf("\n\n")
    
    %check for convergence
    if resnorm(k) < 10e-5
        fprintf("Solution converged!!!\n")
        break
    end

    %check for divergence
    if resnorm(k) > 10e+15
        fprintf("Solution diverged!!!\n")
        break
    end

end

%print solution
fprintf("Spring length vector solution\n")
fprintf("-------------------------------\n")
fprintf("%.5f\n", x(:,1)) %insert solution vector

fprintf("\n\nResidual Norm at every iteration\n")
fprintf("-------------------------------\n")
fprintf("%12.5e \n", transpose(resnorm)) %transpose horizontal vector (resnorm)

function [J] = Jacobian(i, j, x)
    
    delta = 10e-5; %finite diff step size
    
    backx = x; %setting vector equal to guess vector
    backx(j,1) = backx(j,1) - delta; %change one element in the vector

    J = (springfuncs(i, x) - springfuncs(i, backx)) / (delta); %formula for backward finite diff w/ respect to x(j)
    %call f values using springfuncs function


end

function [f] = springfuncs(i, x)

    %define constants
    L = 1; %m
    W = 0.435; %m
    LT = 15; %m
    ki = [1 2 3 4]; %N/m
    kii = [0.125 0.255 0.325 0.425]; %N/m^3 
    % pay clost attention to kii

    %for linear case
    %kii = zeros(4, 1); %N/m^3

    if i == 1
        f  = (x(1)+L) + (x(2)+L) + (x(3)+L) + (x(4)+L) + 3*W - LT;
    end

    if i == 2
        f = ki(1)*x(1) + kii(1)*((x(1))^3) - (ki(2)*x(2)+kii(2)*(x(2))^3);
    end

    if i == 3
        f = ki(2)*x(2) + kii(2)*((x(2))^3) - (ki(3)*x(3)+kii(3)*(x(3))^3);
    end

    if i == 4
        f = ki(3)*x(3) + kii(3)*((x(3))^3) - (ki(4)*x(4)+kii(4)*(x(4))^3);
    end
    %etc
end
