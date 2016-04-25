function [ A L U P time error ] = BLAS3LUPP(A,block_size)
% Block LU factorization with partial pivoting, overwriting L and U on A
% see ALGORITHM 2.10 in Applied Numerical Linear Algebra, J. Demmel, SIAM
% (2007), p.74. Size of blocks is b.
% Uses BLAS2LUPP to produce the LU dec. of rectangular submatrices of A
%
% Altered on 25th of April 2016
% Copyright University of Minho (2016), Computer Science Dpt. 
% HPC Group, Filipe Oliveira and Sergio Caldas
%



n=length(A);

%create diag matrix P
P=eye(n); 
Ainit = A;


start_time = tic;

for i=1:block_size:n-1
    
    [~,p] = max(abs(A(i:n,i)));
    p = p+i-1;

     if p~=i
        A([i p],:) = A([p i], :);
        P( [i p] , : ) =   P( [p i] , : );
     end
    
    % apply row permutations to A and L
    last=min(i+block_size-1,n);
    
    A(i:n,i:last) = BLAS2LUPP(A(i:n,i:last));

    % IF SIZE OF REMAINING BLOCK LARGER THAN b
    if n-i+1 > block_size  
        L22=tril(A(i:last,i:last),-1)+eye(block_size);
        A(i:last,i+block_size:n)=inv(L22)*A(i:last,i+block_size:n); % step 2 (U22)
        A(i+block_size:n,i+block_size:n)=A(i+block_size:n,i+block_size:n)-A(i+block_size:n,i:last)*A(i:last,i+block_size:n); % step 3 (U33)
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
    error = norm ( P*Ainit - L* U ); 

end


