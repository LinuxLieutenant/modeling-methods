function [J] = Jacobian(i, j, x)
    
    delta = 10e-5; %finite diff step size
    
    backx = x; %setting vector equal to guess vector
    backx(j,1) = backx(j,1) - delta; %change one element in the vector

    J = (springfuncs(i, x) - springfuncs(i, backx)) / (delta); %formula for backward finite diff w/ respect to x(j)
    %call f values using springfuncs function


end