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
## @deftypefn {Function File} {@var{retval} =} qr_giv (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: jp <jp@jpordi>
## Created: 2015-04-17

function [Q, R] = qr_giv (R0, decal, rel)
  global epsi
  % xrel contient les relations verticales
  [Dx Dy] = size (R0);
  [Dx nbr] = size (rel);
  q = eye(Dx);
  r = [zeros(Dx, nbr) R0(:, decal+1: Dy)];
  u = rel;
  v = eye(nbr+Dy-decal, nbr);
  [Q, R] = qrupdate (q, r, u, v);
%  spy(abs(Q'*R0)>epsi); grid on; pause
  R = [zeros(Dx, decal) R(:, nbr+1:nbr+Dy-decal)];
 endfunction
