function [ A L U time error ] =BLAS3LU(A,b)
% Block LU factorization with partial pivoting, overwriting L and U on A
% see ALGORITHM 2.10 in Applied Numerical Linear Algebra, J. Demmel, SIAM
% (2007), p.74. Size of blocks is b.
% Uses BLAS2LU to produce the LU dec. of rectangular submatrices of A


start_time = tic;
n=length(A);
Ainit = A;


for i=1:b:n-1
    % apply row permutations to A and L
    last=min(i+b-1,n);  
    
    %Compute the LU factorization with partial pivoting
    A(i:n,i:last)=BLAS2LU(A(i:n,i:last)); % step 1 (L22, L32)
    
    % IF SIZE OF REMAINING BLOCK LARGER THAN b
    if n-i+1 > b  
        L22=tril(A(i:last,i:last),-1)+eye(b);
        A(i:last,i+b:n)=inv(L22)*A(i:last,i+b:n); % step 2 (U22)
        A(i+b:n,i+b:n)=A(i+b:n,i+b:n)-A(i+b:n,i:last)*A(i:last,i+b:n); % step 3 (U33)
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
    error = norm ( Ainit - L* U ) / norm ( Ainit ); 

end


