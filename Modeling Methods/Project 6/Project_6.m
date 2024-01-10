% project 6 least squares
% Project 6 Least Squares 10/27/23

clc, clear, close all;
fprintf("============================================================\n")
fprintf("Project 6\n")
fprintf("\n")
display(datetime())
fprintf("============================================================\n\n")

% Read in the x- y-data
x = importdata("Project 6/X_DATA_FA2023.DAT");
y = importdata("Project 6/Y_DATA_FA2023.DAT");

N = length(x); %length of data vectors

% plot the 
%figure
%scatter(x,y);
%hold on

for order = 1:4 %start:stop 
    if order == 1 || order == 2 || order == 4 %only execute the  code below if order matches the desired polynomial
        figure
        scatter(x,y);
        hold on

        % loop through the order of the least squares fit
        M = order + 1; % number of polynomial terms
        fprintf("\n==========================================\n")
        fprintf('Number of data points = %5i\n', N)
        fprintf('Number of least-squares coefficients = %5i\n', M)
        fprintf('Order of least-squares polynomial fit = %5i\n\n', order)
        
        [coeffs] = LeastSquares(x, y, M);
        fprintf('Least Squares fit coefficients = \n')
        disp(transpose(coeffs));
        fprintf('\n')

        %make x vector for graph
        xrange = 0:0.01:2; %start:step:stop
        xrange = (xrange)';
        ygraph = evalFit(coeffs, xrange); %y values from xgraph

        % evaluate the least squares fit at an array of x values
        %y_lss = evalFit(coeffs, xgraph);
        y_ls = evalFit(coeffs, x); %find the polynomial values of x vector
        
        % sum of squares of errors
        SSE = sum((y-y_ls).^2);
        
        % standard error
        std_error = sqrt(SSE/(N-M));
        
        % R squared
        y_mean = mean(y); %use MATLAB function
        So = sum((y-mean(y)).^2); %aka SSD
        r2 = (So-SSE)/So;
        
        fprintf('Sum of Squares of Errors (SSE) = %12.5e \n', SSE)
        fprintf('Mean= %12.5e \n', y_mean)
        fprintf('R square = %.6f \n', r2)
        fprintf('Standard error= %12.5e \n\n', std_error)
        
        % evaluate the fit at a single x
        x_eval = 1.62;
        y_eval = evalFit(coeffs, x_eval); %polynomial value at x_eval
        
        fprintf('\nEvaluating the LS fit at x = %.6e\n', x_eval)
        fprintf('  SE = %.6e\n',std_error)
        fprintf('   y = %.6e\n',y_eval)
        fprintf('y+SE = %.6e\n',y_eval+std_error) %34% in each direction
        fprintf('y-SE = %.6e\n',y_eval-std_error) %aka 68% together
        fprintf('y+2SE = %.6e\n',y_eval+2*std_error) %47% in each direction 
        fprintf('y-2SE = %.6e\n',y_eval-2*std_error) %aka 95% total

        % plot LS fit and orignal data and LS point on the same plot
        plot(x, y_ls) %polynomial
        plot(x, y_ls+std_error) %+SE
        plot(x, y_ls-std_error) %-SE
        scatter(x_eval, y_eval, 'k', 'filled') %plot x_eval y_eval
        scatter(x_eval, y_eval+std_error, 'k', 'filled', '^') %y_eval +SE
        scatter(x_eval, y_eval-std_error, 'k', 'filled', '^') %y_eval -SE
        grid on
        title('Plot of LS fit of data')
        xlabel('x')
        ylabel('LS fit value y(x)')
        legend("Data Points", "y_L_S(x)", "y_L_S(x)+SE", "y_L_S(x)-SE")
    end
end

function [coeffs] = LeastSquares(x, y, M) % see main for input values

    %Build least squares matrix and RHS vector
    C = zeros(M,M);
    b = zeros(M,1);

    for i = 1:M 
        for j = 1:M 
            C(i,j) = sum(x.^(i+j-2));
            b(i) = sum(y.*x.^(i-1)); %make sure to do element wise .* .^
        end
    end

    %Solve for coeffs
    coeffs = transpose(C\b);
    

end

function [y] = evalFit(a, x)
%import coefficients and x values and return least squares polynomial values
%coefficient vector is called "a" in this script

    %harder but more efficient way:
    i = 1:(length(a)); %vector og exponential powers
    y = sum(a.*x.^(i-1), 2);
    %put x vector to the power of i vector to make a matrix
    %element wise multiply x matrix by horizontal coeffs vector
    %sum each row in the new matrix using a type 2 sum
    %sum ( ,2)


    %easier way: 
    %y = zeros(length(x), 1);

    %for i = 1:length(x) %loop through elements of y
        %for j = 1:length(a) %exponential powers vector

            %y(i) = y(i) + a(j).*x(i).^(j-1); %increase y(i) by current term of a*x
        %end
    %end

end
