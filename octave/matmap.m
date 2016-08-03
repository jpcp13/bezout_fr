function fX = matmap(f, X)
	global n dim
	Id = eye(dim);
	degrees = f{1};
	coeffs = f{2};
	nm = length(coeffs);
	fX = zeros(dim);
	for k = 1:nm
		P = coeffs(k)*Id;
		for j = 1:n
			P = P*X(:, :, j)^degrees(k, j);
		endfor
		fX = fX + P;
	endfor
endfunction
