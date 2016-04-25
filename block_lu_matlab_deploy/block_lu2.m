%The MIT License (MIT)

%Copyright (c) 2014 sgh1.net

%Permission is hereby granted, free of charge, to any person obtaining a copy
%of this software and associated documentation files (the "Software"), to deal
%in the Software without restriction, including without limitation the rights
%to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%copies of the Software, and to permit persons to whom the Software is
%furnished to do so, subject to the following conditions:

%The above copyright notice and this permission notice shall be included in
%all copies or substantial portions of the Software.

%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
%THE SOFTWARE.


function [LU,piv] = block_lu2(A,blocksize)

	n = length(A);
	
	%setup permuation vector
	piv = zeros(1,n);
	for kk = 1:n
		piv(kk) = kk;
	end

	for k = 1:blocksize:n
		
		vb = min([blocksize,n - k + 1])
		
		%get cur block of A, and do normal LU on it
		[piv,LU] = DoolittleLU_block2(A,(n-k+1),((k+vb-1)-k)+1,k,k,k-1,piv);
		A = LU;

		%do inter-block update, some expensive updates
		%here, but remember it's just one block
		A_sub = A(k:k+vb-1, k:k+vb-1);
        
		AA = tril_( A_sub ) + eye(vb);
		temp3 = AA^-1 * A(k:k+vb-1, k+vb:n);
		A(k:k+vb-1, k+vb:n) = temp3;

		%do submatrix mul update
		A(k+vb:n,k+vb:n) = A(k+vb:n,k+vb:n) - A(k+vb:n,k:k+vb-1) * A(k:k+vb-1,k+vb:n);
		

	end

	disp('lu decomposition done');
	LU = A

end
