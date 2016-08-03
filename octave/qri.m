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
## @deftypefn {Function File} {@var{retval} =} qri (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: jp <jp@jpordi>
## Created: 2015-04-12

function [Q, R, P] = qri (R0, decal)
global epsi
  [sx, sy] = size (R0);
  Q = eye(sx); R = R0; P = eye(sy);
  if sx+decal < sy
    idy = [sx+decal+1: sy, 1: sx+decal];
    R = R(:, idy); P = P(idy, :);
    Qt = P'; Rt = R'; Pt = Q';
    Qt = fliplr(Qt); Rt = flipud(Rt); 
    Rt = fliplr(Rt); Pt = flipud(Pt);
    [q, r, p] = qr(Rt(sx+1: sy, :));
%    nbr = sum (abs (diag (r)) > epsi);
    nbr = rank(r);
    Qt(:, sx+1:sy) = Qt(:, sx+1:sy)*q;
    Rt(sx+1: sy, :) = q'*Rt(sx+1: sy, :);
    qt = eye(sy);
    rt = Rt; rt(sx+1: sy, :) = 0;
    u = zeros(sy, nbr); u(sx+1: sx+nbr, :) = eye(nbr);
    v = Rt(sx+1: sx+nbr, :)';
    [q1, Rt] = qrupdate (qt, rt, u, v);
    Qt = Qt*q1;
    Rt = fliplr(Rt); Pt = flipud(Pt);
    Qt = fliplr(Qt); Rt = flipud(Rt);
    Q = Pt'; R = Rt'; P = Qt';
  endif
endfunction
