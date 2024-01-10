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