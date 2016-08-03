function id = lin2mult(u, deg)
n = length(deg);
p = prod(deg);
id = [];
for j = 1:n
	p = p/deg(j);
	q = floor(u/p);
	id = [id q];
	u = u - q*p;
endfor
