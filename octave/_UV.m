function [U, V] = _UV()
global n deg
U = cell(1, n);
V = cell(1, n);
for j = 1:n
	U{j} = fft([0 1 zeros(1, j*deg(j)-2)]);
	V{j} = fft([0 1 zeros(1, (n-j+1)*deg(j)-2)]);
endfor
