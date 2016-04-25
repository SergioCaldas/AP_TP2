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


function [piv_out,LU] = DoolittleLU_block2(A,nrows,ncols,row_offset,col_offset,piv_offset,piv_in)
	
	colend = ncols + col_offset -1
	rowend = nrows + row_offset -1			
				
	for k = col_offset:colend

		%find pivot row
		piv_val = 	abs( A(k,k) );
		piv_in(k) = k;
		
		for kk = k:nrows
			if(	abs( A(kk,k) ) > piv_val )
				piv_val = 	abs( A(kk,k) );
				piv_in(k) = kk;
			end
		end
				
		%do complete row swap -- this will reach outside
		%of the block currently being worked on!
		temp = A( piv_in(k) , : );
		A(piv_in(k), : ) = A( k , : );
		A( k, : ) = temp; 	
	
		A( k+1:rowend, k ) = A( k+1:rowend, k )./ A(k,k); 
				
		%do the gauss-like operation
		for ii = k+1:rowend
			for jj = k+1:colend
				A(ii,jj) = A(ii,jj) - A(k,jj)*A(ii,k); 
			end
		end			
	end

	LU = A;	
	piv_out = piv_in;
	
	
