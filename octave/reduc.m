## Copyright (C) 2015 jp
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {Function File} {@var{retval} =} reduc (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: jp <jp@jpordi>
## Created: 2015-04-25

function [B] = reduc (B, decal, j)
	global epsi n
	[Dx, Dy] = size(B(:, :, 1));
	idx = 1:Dx; 
	idy = 1:decal;
	rel_j = B(idx, idy, j+1);
	[q, r, p] = qr(rel_j');
	nbr = sum (abs (diag (r)) > epsi); % nbr = nombre de nouvelles relations
	B = apply (eye(Dx), q', B, idx, idy);
	rel_j = B(:, 1:nbr, j+1);
	[Q R] = qr_giv (B(:, :, 1), decal, rel_j);
	B = apply (Q, eye(Dy), B, 1:Dx, 1:Dy);
	B(1:nbr, :, :) = [];
endfunction
