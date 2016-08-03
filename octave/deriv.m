function dfj = deriv(f, j)
	degrees = f{1};
	coeffs = f{2};
	degj = degrees(:, j);
	coeffs = coeffs .* degj';
	degrees(:, j) = max(0, degj - 1);
	dfj = {degrees, coeffs};
endfunction

