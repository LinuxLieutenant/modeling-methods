function [f] = springfuncs(i, x)

    %define constants
    L = 1; %m
    W = 0.435; %m
    LT = 15; %m
    ki = [1 2 3 4]; %N/m
    %kii = [0.125 0.255 0.325 0.425]; %N/m^3 
    % pay clost attention to kii

    %for linear case
    kii = zeros(4, 1); %N/m^3

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