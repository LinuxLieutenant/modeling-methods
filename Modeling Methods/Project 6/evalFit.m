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