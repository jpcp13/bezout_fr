function y = map(f, x)
	degrees = f{1};
	coeffs = f{2};
	nm = length(coeffs);
	ox = ones(nm, 1)*x;
	p = prod(ox.^degrees, 2);
	y = coeffs*p;
endfunction
