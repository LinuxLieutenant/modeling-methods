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
