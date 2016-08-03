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
## @deftypefn {Function File} {@var{retval} =} fourier (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: jp <jp@jpordi>
## Created: 2015-05-23

function [Fx, Fy] = fourier()
  global n degx degy
  Fx = 1;
  Fy = 1;
  for j = 1:n
    Fxj = fft(eye(degx(j)));
    Fyj = fft(eye(degy(j)));
    Fx = kron(Fx, Fxj);
    Fy = kron(Fy, Fyj);
  endfor
endfunction
