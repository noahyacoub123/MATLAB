%Noah Yacoub, 21026083
function out = naive_gauss(A,b)
%Naive Gaussian elimination given a matrix and the solution vector, without pivoting 

[row, col] = size(A);
n = row;
if row ~= col
    error("Matrix [A] is not a square matrix") %Checks if the matrix is square
end

if row ~= length(b)
    error("Matrix dimension mismatch") %Checks that the dimentions of the solutions 
end
    
if det(A) == 0
    error("Matrix [A] is a singular matrix") %Makes sure that the matrix isn't singular 
end

a = [A b]; %creating a new, augmented matrix of the given matrix and solution vector

for k = 1:row-1
    for j = k+1:row
        a(j,k:n+1) = a(j,k:n+1) - a(j,k)/a(k,k)*a(k,k:n+1);   %The forward elimination part of the process 
    end
end

x = zeros(n, 1);
x(n) = a(n, n+1) / a(n, n);


for k = row-1:-1:1
    
       x(k) = (a(k, n + 1) - a(k, k + 1:n)*x(k+1:n)) / a(k, k);  %The backwards substitution part of the process 

end

out = x;
end
