% LU factorization with partial pivoting, overwriting L and U on A
% see ALGORITHM 2.9 in Applied Numerical Linear Algebra, J. Demmel, SIAM
% (2007), p.72
% 
% Altered on 25th of April 2016
% Copyright University of Minho (2016), Computer Science Dpt. 
% HPC Group, Filipe Oliveira and Sergio Caldas
%

function [ A P ]  = BLAS2LUPP(A ,  P )

[m n] = size(A);
L = eye(n);
U = A;
Ainit = A;

%create diag matrix P

start_time = tic;
for i=1:min(m-1,n)

    [~,p] = max(abs(A(i:n,i)));
    p = p+i-1;

    % apply row permutations to A and P
    if p~=i
        A([i p],:) = A([p i], :);
        P( [i p] , : ) =   P( [p i] , : );
    end
    A(i+1:m,i)=A(i+1:m,i)/A(i,i);
    if i<n
        A(i+1:m,i+1:n)=A(i+1:m,i+1:n)-A(i+1:m,i)*A(i,i+1:n);
    end

end
duration = toc(start_time);


