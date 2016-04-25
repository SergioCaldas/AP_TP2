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

function [x] = LU_Backsolve(LU,b,pivot)

	%backsolve code, nothing special here
	%just  basic backsolve of an LU
	%decomposition with partial pivoting

	%there may be some room for optimization
	
	n = length(LU);
	x = zeros(1,n) ;
	
	for k = 1:n
		if( pivot(k) ~= k )
			temp = b(k);
			b(k) = b( pivot(k) );
			b( pivot(k) ) = temp;
		end
		
		x(k) = b(k);
		
		for ii = 1:k-1
			%ii
			x(k) = x(k) - x(ii)*LU(k,ii); 
		end
	end

	
	for k = n:-1:1
		if( pivot(k) ~= k )
			temp = b(k);
			b(k) = b( pivot(k) );
			b( pivot(k) ) = temp;
		end

		for ii = k+1:n
			x(k) = x(k) - x(ii)*LU(k,ii);
		end
		
		x(k) = x(k) / LU(k,k);
		
	end

end
