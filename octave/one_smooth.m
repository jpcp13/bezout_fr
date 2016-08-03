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
## @deftypefn {Function File} {@var{retval} =} smooth (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: jp <jp@jpordi>
## Created: 2015-04-26

function [B, decal] = one_smooth (B, decal)
  global epsi
  [Dx, Dy] = size(B(:, :, 1));
  d = diag (B(:, :, 1), decal);
  ld = length(d);
  z = find (abs(d) < epsi);
  idz = min (z);
  nz = length (z);
  if nz > 0 
    % 1er bloc
    idx = 1: idz - 1;
    idy = 1: decal + idz;
    r0 = B(idx, idy, 1);
    [q, r, p] = qri (r0, decal);
    B = apply(q, p, B, idx, idy);
    % 2eme bloc
    idx = idz: -decal + Dy;
    idy = decal + idz + 1: Dy;
    q0 = eye(-decal + Dy - idz);
    r0 = B(idx, idy, 1);
    p = eye(-decal + Dy - idz);
    [q, r] = qrinsert (q0, r0(2: end, :), 1, r0(1, :), 'row');
    B = apply(q, p, B, idx, idy);
    decal = decal + 1;
    endif
  endfunction
