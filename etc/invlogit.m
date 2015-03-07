function [ x ] = invlogit(x)
%INVLOGIT	Calculate the inverse logit of the input
%	[ x ] = invlogit(x) calculates the inverse logit of the input:
%       x = 1./(1+exp(-x))
%   Usually it is used to convert a value into a probability.
%   
%	Inputs:
%		x   - Vector or matrix of values between -Inf and Inf
%		
%
%	Outputs:
%		x   - Vector or matrix of values between 0 and 1
%		
%
%	Example
%       x = round(-50 + 100*rand(100,1));
%		[ x ] = invlogit(x); 
%	
%	See also LOGIT

%	References:
%	
%	

%	Copyright 2012 Alistair Johnson
%	This program is free software: you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation, either version 3 of the License, or
%	(at your option) any later version.
%	
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%	
%	You should have received a copy of the GNU General Public License
%	along with this program.  If not, see <http://www.gnu.org/licenses/>.


%	$LastChangedBy$
%	$LastChangedDate$
%	$Revision$
%	Originally written on GLNXA64 by Alistair Johnson, 10-Aug-2012 13:53:14
%	Contact: alistairewj@gmail.com

x = 1./(1+exp(-x));

end