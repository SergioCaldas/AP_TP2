A = rand(6); 
%[L,U,P] = lu(A);

%n=norm(A,'inf');

A = 10*A;
[L,U,P] = lu(A);

C = P*A-L*U;