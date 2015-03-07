function [ x ] = logit(x)
%LOGIT	Calculate the logit of the input
%	[ x ] = logit(x) calculates the logit of the input:
%       x = log(x./(1-x))
%   Usually it is used to convert a probability into a value.
%
%	Inputs:
%		x   - A vector or matrix of values between 0 and 1
%
%	Outputs:
%		x   - A vector or matrix of values between -Inf and Inf
%		
%
%	Example
%       x = rand(100,1);
%		[ x ] = logit(x);
%	
%	See also INVLOGIT

%	Copyright 2012 Alistair Johnson

%	$LastChangedBy$
%	$LastChangedDate$
%	$Revision$
%	Originally written on GLNXA64 by Alistair Johnson, 10-Aug-2012 13:53:05
%	Contact: alistairewj@gmail.com
x(x==0) = x(x==0)+1e-10;
x(x==1) = x(x==1)-1e-10;
x = log(x./(1-x));

end
