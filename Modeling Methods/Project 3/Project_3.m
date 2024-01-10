%
%Project 3 Gauss Elimination
%EML 3034 Modeling Methods
%Due 9-29-23

clc, clear, close all

A= dlmread("A-6.txt"); %make sure file name matches
b= dlmread("b-6.txt"); %File must be in working directory

format shortE

%% Call diag dom
diag_dom(A);

%% Solution 1: Gauss Elimination
[x,Atri,btri] = gauss(A,b);


%print A and b after triangulation
fprintf("A(10,50) = %.4e\n", Atri(10,50)) %Change the range to whatever the quiz says (4,10)
fprintf("b(85) = %.4e\n", btri(85)) %also change this to what quiz says

%print solution elements
fprintf("\nGauss Answers\n")    %also change thses values to what the quiz says
fprintf("x(5) = %.4e\n", x(5))
fprintf("x(10) = %.4e\n", x(10))
fprintf("x(25) = %.4e\n", x(25))
fprintf("x(75) = %.4e\n", x(75))

%Calcultate the residual
%infinity norm
res = A*x - b;
res_norm = norm(res, Inf)

%% Solution 2 MATLAB
x_m = linsolve(A,b);
%print solution elements
fprintf("\nGauss Answers\n")    %also change thses values to what the quiz says
fprintf("x_m(5) = %.4e\n", x_m(5))
fprintf("x_m(10) = %.4e\n", x_m(10))
fprintf("x_m(25) = %.4e\n", x_m(25))
fprintf("x_m(75) = %.4e\n", x_m(75))

%For quiz (q4):
A= dlmread("A-6.txt"); %make sure file name matches
b= dlmread("b-6.txt"); %File must be in working directory
%
x = ones(length(b),1);
%
res = A*x - b;
%
quiz_norm = norm(res, Inf)
%
%

function [x, A, b] = gauss(A,b)
    fprintf("Converting A to upper triangular:\n")
    n = length(A); % get number of equations
    x = zeros(length(b),1); % allocate memory

    % Triangulate the A matrix
    for j=1:n-1 % ?loop through each column that has elements below diagonal
        for i=j+1:n % ?transform all elements below diagonal to 0
            s = -A(i,j)/A(j,j); % calculate multiplier
            A(i,:) = A(i,:) + s*A(j, :); % apply ERO to the row of A
            b(i) = b(i) + s*b(j); % apply ERO to b (right hand side vector) as well
        end
    end

    % solve last equation
    x(end) = b(end) / A(end, end);

    % solve subsequent equations moving backwards (up)
    for i= (n-1) : -1 : 1
        x(i) = (b(i) - A(i, (i+1) : n) * x((i+1) : n)) / A(i, i);
    end

end

function diag_dom(A)

    n = length(A);
    dom = 1;

    for i = 1:n
        if ( sum(abs(A(i, :))) - A(i,i) ) < A(i,i) %diagonal is less than the rest of the row
            dom = 0;
            fprintf("\nMatrix is NOT strictly diagonally dominant\n")
            break
        end
    end
    if dom == 1
            fprintf("\nMatrix IS strictly diagonally dominant\n")
    end
end
