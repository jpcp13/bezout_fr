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
## @deftypefn {Function File} {@var{retval} =} triang_block_giv (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: jp <jp@jpordi>
## Created: 2015-04-25

function [B_out, decal] = triang_block_giv (B, bls)
	global n
	blsx0 = bls;
	blsy0 = bls;
	R0 = B(:, :, 1);
	[Q, R, P, blsx, blsy] = qr_block_giv (R0, blsx0, blsy0);
	norm_test = norm(Q*R*P - R0)
	decal = sum(blsy) - sum(blsx);
	B_out = zeros(size(B));
	for j = 0:n
		B_out(:, :, j+1) = Q'*B(:, :, j+1)*P';
	endfor
 endfunction
