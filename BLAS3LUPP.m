function [ A time error ] = BLAS3LUPP(A, block_size)
% Block LU factorization with partial pivoting, overwriting L and U on A
% see ALGORITHM 2.10 in Applied Numerical Linear Algebra, J. Demmel, SIAM
% (2007), p.74. Size of blocks is b.
% Uses BLAS2LUPP to produce the LU dec. of rectangular submatrices of A
%
% Altered on 25th of April 2016
% Copyright University of Minho (2016), Computer Science Dpt. 
% HPC Group, Filipe Oliveira and Sergio Caldas
%

[m n] = size(A);

%create diag matrix P
 
Ainit = A;

P = eye(n);
start_time = tic;

for i=1:block_size:n-1

    for row_i=i:i+block_size-1
    [~,p] = max(abs(A(row_i:m,row_i)));
    
    p=p+row_i-1;
    % apply row permutations to A and P
    if p~=row_i
        A([row_i p],:) = A([p row_i], :);
        P([row_i p],:) = P([p row_i], :);
    end
    end
    
    last=min(i+block_size-1,n);  
    
    % Apply LU with pivoting to A(k:n,k:k+vb?1)    
    [A__new , Lu, Pu,  Pi, time, error ] = BLAS2LUPP( A(i:m,i:last));
   
    A(i:m,1:n) = Pi * A(i:m,1:n) ;
    P(i:m,1:n) = Pi * P(i:m,1:n) ;
    A(i:m,i:last) = A__new;
    

    % IF SIZE OF REMAINING BLOCK LARGER THAN b
    if n-i+1 > block_size  
        L22=tril(A(i:last,i:last),-1)+eye(block_size);
        A(i:last,i+block_size:n)=inv(L22)*A(i:last,i+block_size:n); % step 2 (U22)
        A(i+block_size:n,i+block_size:n)=A(i+block_size:n,i+block_size:n)-A(i+block_size:n,i:last)*A(i:last,i+block_size:n); % step 3 (U33)
    end
end

duration = toc(start_time);

if nargout > 2
   time = duration;
   L=tril(A);
   U=triu(A);
   
   [ linesL colsL ] = size ( U );

    for pos=1:linesL
    L ( pos,pos ) = 1;
    end
    
    error = norm ( P * Ainit - L * U ) / norm(Ainit); 

end

