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
