function [ A L U time error ] =BLAS2LU(A)
% LU factorization with partial pivoting, overwriting L and U on A
% see ALGORITHM 2.9 in Applied Numerical Linear Algebra, J. Demmel, SIAM
% (2007), p.72
% This differs from LUfact1 because A may not be square, has m rows and n
% columns, with m>=n
%


[m n]=size(A);
start_time = tic;
Ainit = A;

for i=1:min(m-1,n)
    % apply row permutations to A and L
    A(i+1:m,i)=A(i+1:m,i)/A(i,i);
    if i<n
        A(i+1:m,i+1:n)=A(i+1:m,i+1:n)-A(i+1:m,i)*A(i,i+1:n);
    end
end
    
duration = toc(start_time);
if nargout > 1
   time = duration;
    L=tril(A);
   U=triu(A);
   [ linesL colsL ] = size ( U );

    for pos=1:linesL
    L ( pos,pos ) = 1;
    end
    Ainit = A;
   error = norm ( Ainit - L* U ) / norm ( Ainit ); 
   
end

