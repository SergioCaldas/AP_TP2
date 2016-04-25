function [ L U P time ] = BLAS3LUPP(A,block_size)
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
L = eye(n);
U = A;
%create diag matrix P
P=eye(n); 
start_time = tic;

for i=1:block_size:n-1
    
    % apply row permutations to A and L
    last=min(i+block_size-1,n);
    
    [ L_block U_block P_block ] = BLAS2LUPP(A(i:n,i:last));
   % A(i:n,i:last) = L_block * U_block;

    % SIZE OF REMAINING BLOCK LARGER THAN b
    if n-i+1 > block_size  
        L22=L_block+eye(block_size);
        L(i:last,i+block_size:n)=inv(L22)*U(i:last,i+block_size:n); % step 2 (U22)
        U(i+block_size:n,i+block_size:n)=U(i+block_size:n,i+block_size:n)-L(i+block_size:n,i:last)*U(i:last,i+block_size:n); % step 3 (U33)
    end
end

duration = toc(start_time);

if nargout > 1
   time = duration;
end

