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
## @deftypefn {Function File} {@var{retval} =} qr_block (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: jp <jp@jpordi>
## Created: 2015-04-14

function [Q, R, P, blsx, blsy] = qr_block_giv (R0, blsx0, blsy0)
  global epsi
  Blsx0 = cumsum(blsx0);
  [Dx, Dy] = size (R0);
  Q = eye(Dx); R = R0; P = eye(Dy);
  sx = 0; sy = 0; blsx = 0*blsx0;
  decal = 0;
  for j = 1: length (blsx0)
    idx = sx+1: Blsx0(j);
    idy = sy+1: sy+blsy0(j);
    [q, r, p] = qr(R(idx, idy)); % rank revealing QR !!
    ds = sum(abs(diag(r)) > epsi);
    Q(:, idx) = Q(:, idx)*q;
    R(idx, :) = q'*R(idx, :);
    R(:, idy) = R(:, idy)*p;
    P(idy, :) = p'*P(idy, :);
    sx = sx + ds;
    sy = sy + blsy0(j);
    if sx > 0
      [q, r, p] = qri (R(1:sx, 1:sy), decal);
      Q(:, 1:sx) = Q(:, 1:sx)*q;
      R(1:sx, :) = q'*R(1:sx, :);
      R(:, 1:sy) = R(:, 1:sy)*p';
      P(1:sy, :) = p*P(1:sy, :);
    endif
    decal = decal + blsy0(j) - ds;
    blsx(j) = ds;
  endfor
  blsy = blsy0;
endfunction
