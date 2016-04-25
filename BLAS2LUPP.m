% LU factorization with partial pivoting, overwriting L and U on A
% see ALGORITHM 2.9 in Applied Numerical Linear Algebra, J. Demmel, SIAM
% (2007), p.72
% 
% Altered on 25th of April 2016
% Copyright University of Minho (2016), Computer Science Dpt. 
% HPC Group, Filipe Oliveira and Sergio Caldas
%

function [ L U P ]  = BLAS2LUPP(A)

[m n] = size(A);

L = eye(n);
U = A;
%create diag matrix P
P=eye(n); 

for i=1:min(m-1,n)

    [~,p] = max(abs(A(i:n,i)));
    p = p+i-1;

    % apply row permutations to A and P
    if p~=i
        U([i p],:) = U([p i], :);
        P( [i p] , : ) =   P( [p i] , : );
        L([i,p],1:p-1) =  L([i,p], 1:p-1);
    end

   for h = i+1:n      
     L(h, i) = U(h, i) / U(i, i);
     U(h, :) =  U(h, :) - L(h, i)*U(i, :);
   end

end
 
end
