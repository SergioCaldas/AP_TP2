function A =BLAS2LUPP(A)
% LU factorization with partial pivoting, overwriting L and U on A
% see ALGORITHM 2.9 in Applied Numerical Linear Algebra, J. Demmel, SIAM
% (2007), p.72
% This differs from LUfact1 because A may not be square, has m rows and n
% columns, with m>=n

[m n]=size(A);
L=eye(n); 
P=eye(n); 
U=A;
for i=1:min(m-1,n)
    [~,p] = max(abs(A(i:n,i)));
    p = p+i-1;
    if p~=i
        A([i p],:) = A([p i], :);
    end
    % apply row permutations to A and L
    A(i+1:m,i)=A(i+1:m,i)/A(i,i);
    if i<n
        A(i+1:m,i+1:n)=A(i+1:m,i+1:n)-A(i+1:m,i)*A(i,i+1:n);
    end
end
end
