%a = reshape(1:24, [2 3 4])
%save -6 octave_a.mat a

load np_B.mat

global epsi n dim
epsi = 1e-7
[Dx, Dy, n1] = size (B);
n = n1 - 1;
for j = 0:n
	B(:, :, j+1) = fliplr(flipud(B(:, :, j+1)'));
endfor
bls = fliplr(bls);
%degrees, coeffs, nb_monomials
f = cell(1, n);
first = 1 + cumsum([0 nb_monomials(1:n-1)]);
last = cumsum(nb_monomials);
for i = 1:n
	f{i} = {degrees(first(i):last(i), :), coeffs(first(i):last(i))};
endfor
df = cell(n, n);
for i = 1:n
	for j = 1:n
		df{i, j} = deriv(f{i}, j);
	endfor
endfor

B0 = B;


[B, decal] = triang_block_giv (B, bls);
%plot(log10(eps + abs(diag(B(:, :, 1), decal) ) ), 'r*-'); grid; pause

delta = 1;
while delta > 0
  [B, decal1] = grand_reduc (B, decal);
  delta = decal1 - decal;
  decal = decal1
endwhile
dim = Dy - decal

%%%%%%

%B(dim+1:end, :, :) = [];
%B(:, 1:decal, :) = [];
%B(:, :, 1) = triu(B(:, :, 1));

%X = zeros(dim, dim, n);
%Xred = zeros(dim, dim, n);
%for j = 1:n
%	X(:, :, j) = B(:, :, 1)\B(:, :, j+1);
%	Xred(:, :, j) = Bred(:, :, 1)\Bred(:, :, j+1);
%endfor

%Bchow = zeros(dim, dim);
%Bredchow = zeros(dim, dim);
%for j = 1:n
%	Bchow = Bchow + rand*B(:, :, j+1);
%	chow = B(:, :, 1)\Bchow;
%	Bredchow = Bredchow + rand*Bred(:, :, j+1);
%	redchow = Bred(:, :, 1)\Bredchow;
%endfor

%[W, S] = schur(chow, "complex");
%[U, T] = schur(redchow, "complex");
%rac_schur = zeros(dim, n);
%rac_red = zeros(dim, n);
%for j = 1:n
%	rac_schur(:, j) = diag(W'*X(:, :, j)*W);
%	rac_red(:, j) = diag(U'*Xred(:, :, j)*U);
%endfor

%f_rac_schur = zeros(n, dim);
%f_rac_red = zeros(n, dim);
%for i = 1:n
%	for k = 1:dim
%		f_rac_schur(i, k) = map(f{i}, rac_schur(k, :));
%		f_rac_red(i, k) = map(f{i}, rac_red(k, :));
%	endfor
%endfor
%subplot(2, 1, 1)
%plot(log10(eps + abs(f_rac_schur')), '*-'); grid
%subplot(2, 1, 2)
%plot(log10(eps + abs(f_rac_red')), '*-'); grid


