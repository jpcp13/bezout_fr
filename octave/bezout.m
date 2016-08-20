load np_B.mat

global n epsi
n = size(B, 3) - 1
D = size(B, 1)
epsi = 1e-9
%
%%deg = [2 3 4]
%n = length(deg)
%%nb_monomials = 15*ones(1, n);
%nbm = sum(nb_monomials)
%%degrees = [];
%%for j = 1:n
%%  degrees = [degrees randi([0, deg(j)], nbm, 1)];
%%endfor
%%coeffs = randi([-10 10], 1, nbm);
%
%%degrees = [2 0; 1 2; 0 0; 2 1; 1 0]
%%coeffs = [1 1 -1 1 1]
%%nb_monomials = [3 2]
%nbm = sum(nb_monomials);
%deg = max(degrees);
%n = length(deg)
%
%degx = (1:n).*deg;
%degy = (n:-1:1).*deg;
%D = prod(degx)
%
%f = cell(1, n);
%first = 1 + cumsum([0 nb_monomials(1:n-1)]);
%last = cumsum(nb_monomials);
%for i = 1:n
%	f{i} = {degrees(first(i):last(i), :), coeffs(first(i):last(i))};
%endfor
%df = cell(n, n);
%for i = 1:n
%	for j = 1:n
%		df{i, j} = deriv(f{i}, j);
%	endfor
%endfor
%
%[U, V] = _UV();
%[Fx, Fy] = fourier();
%
%
%[ax, ay] = _axy();
%%B = zeros(D, D, n+1);
%%for k = 0:n
%%	tic(); C = _mat(k); toc()
%%	Bk = Fx\C/Fy;
%%	Bk = real(Bk);
%%	Bk = round(Bk);
%%	B(:, :, k+1) = fliplr(Bk(ax, ay));
%%endfor
%
%bls = [];
%for j = 1:n
%  bls = [bls D/n/deg(j)*ones(1, deg(j))];
%endfor

B0 = B;
tic
[B, decal] = triang_block_giv (B, bls);
triang_time = toc
save -ascii octave_triang_time.txt triang_time

dbt = diag(B(:,:,1), decal);
save -ascii diag_beztri.txt dbt

tic
delta = 1;
while delta > 0
  [B, decal1] = grand_reduc (B, decal);
  delta = decal1 - decal;
  decal = decal1
endwhile
dim = D - decal
reduct_time = toc
csvwrite ('octave_reduct_time.txt', round(1000*reduct_time))

B(dim+1:end, :, :) = [];
B(:, 1:decal, :) = [];
B(:, :, 1) = triu(B(:, :, 1));
dbtf = diag(B(:,:,1), 0);
save -ascii diag_beztri_final.txt dbtf

X = zeros(dim, dim, n);
for j = 1:n
	X(:, :, j) = B(:, :, 1)\B(:, :, j+1);
endfor


[W, S] = schur(X(:, :, 1), "complex");
rac = zeros(dim, n);

for j = 1:n
	rac(:, j) = diag(W'*X(:, :, j)*W);

endfor

rr = real(rac);
ir = imag(rac);
save -ascii rac_real.txt rr
save -ascii rac_imag.txt ir


%f_rac = zeros(n, dim);
%for i = 1:n
%	for k = 1:dim
%		f_rac(i, k) = map(f{i}, rac(k, :));
%	endfor
%endfor
%
%
%plot(log10(eps + abs(f_rac')), '*-'); grid


