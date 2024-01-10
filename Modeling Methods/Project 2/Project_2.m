% 
% Project 2 Newton Rapshon and Secant Method
% EML 3034H
% Due Date: 9/15/23
clc, clear, close all;
fprintf("============================================================\n")
fprintf("Project 2 - finding the root of a scalar equation\n")
fprintf("\n")
display(datetime("today"))
fprintf("============================================================\n\n")

L = 2500;
D = 1;
r = 3 * 10^-4;
V = 0.35;
Re = 3.1818 * 10^4;
g = 32.2;

f = @(x) (x.^(-0.5)) + 2 .* log10((r./3.7)+(2.51./(Re.*x.^(0.5))));
df = @(x) (-0.5).*(x.^(-1.5)) - ( ((2.51./Re) .* log10(exp(1))) ./ ( (r./3.7)+(2.51./(Re.*x.^0.5)) )) .* (x.^(-1.5));
h = @(f) f .* (L./D) .* ( (V.^2) ./ (2.*g) );
max_iter = 1000;
guess1 = 0.001;
guess2 = 0.01;
res_tol = 1e-5;
conv_tol = 1e-5;
% plot the function
x = 0.001:0.001:0.1;
fval = f(x);
figure
plot(x, fval, x, zeros(length(x)))
title('Function to find the root')
% Run the solvers
fprintf('Newton Raphson\n')
% call the newton_raphson() function here
[nr_x, nr_conv, nr_res, nr_iters] = newton_raphson_func(f, df, guess1, conv_tol, res_tol, max_iter);

fprintf('\n\n=====================================================\n')
fprintf('Secant\n')
% call the secant() function here
[sec_x, sec_conv, sec_res, sec_iters] = secant_func(f, guess1, guess2, conv_tol, res_tol, max_iter);

fprintf('\n==========================================================\n')
fprintf('Answers:\n')
fprintf(' newton = %.6e in %i iterations.\n', nr_x(end), nr_iters);
fprintf(' secant = %.6e in %i iterations.\n', sec_x(end), sec_iters);

% Headloss
nr_head = h(nr_x(end));
sec_head = h(sec_x(end));
fprintf("\nNR Head Loss: %.6f\nSec Head Loss: %.6f\n", nr_head, sec_head);

% plot the residual and convergences
plot_iters_nr = 0:nr_iters; %start:stop 0 to num of iterations nr
figure
ax = plotyy(plot_iters_nr, nr_res, plot_iters_nr, nr_conv);
title('Newton-raphson')
xlabel('iterations')
ylabel(ax(1), 'residual')
ylabel(ax(2), 'convergence')
grid on

plot_iters_s = 0:sec_iters+1;
figure
ax = plotyy(plot_iters_s, sec_res, plot_iters_s, sec_conv);
title('secant')
xlabel('iterations')
ylabel(ax(1), 'residual')
ylabel(ax(2), 'convergence')
grid on

function [x, conv, res, i] = newton_raphson_func(f, df, initial_guess, conv_tol, res_tol, max_iter)
    % initial guess
    x(1) = initial_guess;
    % initial convergence
    conv(1) = 0;
    % initial residual
    res(1) = abs(f(x(1)));
    fprintf('iteration x residual convergence\n');
    fprintf('========= ===== ========= ==========\n');
    fprintf(' %2i %.8e %.4e %.4e\n', 0,x(1),res(1),conv(1));
    % initialize the iteration counter
    i = 1;
    % iterate using our formula
    while i<max_iter
        fval = f(x(i));
        dfval = df(x(i));
        if (dfval==0) || (isnan(dfval))
            error("!!! FAILED. The derivative value is %i. Try another guess.", dfval)
            break;
        end
        x(i+1) = x(i)-(fval./dfval); % apply Newton-Raphson formula
        conv(i+1) = abs(x(i+1)-x(i)); % calculate iterative convergence
        res(i+1) = f(x(i+1)); % calculate the residual
        fprintf('    %i    %.8e    %.4e    %.4e\n', ...
                      i,x(i+1),res(i+1),conv(i+1));
        % check for convergence
        if res(i+1)<res_tol && conv(i+1)<conv_tol
            break;
        end
        i=i+1; % increment our iteration counter
    end
end

function [x, conv, res, i] = secant_func(f, initial_guess1, initial_guess2, conv_tol, res_tol, max_iter)
    % initial guess
    x(1) = initial_guess1;
    x(2) = initial_guess2;
    % initial convergence
    conv(1) = 0;
    conv(2) = abs(x(2)-x(1));
    % initial residual
    res(1) = abs(f(x(1)));
    res(2) = abs(f(x(2)));
    fprintf('iteration x residual convergence\n');
    fprintf('========= ===== ========= ==========\n');
    fprintf(' %2i %.8e %.4e %.4e\n', 0,x(1),res(1),conv(1));
    fprintf(' %2i %.8e %.4e %.4e\n', 0,x(2),res(2),conv(2));
    % initialize the iteration counter
    i = 2;
    % iterate using our formula
    while i<max_iter
        fval = f(x(i));
        %dfval = df(x(i));
        %if (dfval==0) || (isnan(dfval))
        %    error("!!! FAILED. The derivative value is %i. Try another guess.", dfval)
        %    break;
        %end
        x(i+1) = x(i) - ( fval.*(x(i)-x(i-1)))./(fval-f(x(i-1))); % apply Secant formula
        conv(i+1) = abs(x(i+1)-x(i)); % calculate iterative convergence
        res(i+1) = abs(f(x(i+1))); % calculate the residual
        fprintf('     %i     %.8e    %.4e    %.4e\n', ...
                       i-1, x(i+1), res(i+1), conv(i+1));
        % check for convergence
        if res(i+1)<res_tol && conv(i+1)<conv_tol
            break;
        end
        i=i+1; % increment our iteration counter
    end
        i=i-1; %
end
