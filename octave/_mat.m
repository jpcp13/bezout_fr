function C = _mat(k)
  global n D degx degy U V
  C = zeros(D, D);
  for iu = 1:D
    for jv = 1:D
      u = 1 + lin2mult(iu-1, degx);	v = 1 + lin2mult(jv-1, degy);
      uu = zeros(1, n); vv = zeros(1, n);
      for j = 1:n
        uu(j) = U{j}(u(j));
        vv(j) = V{j}(v(j));
      endfor
      m = _pol(uu, vv, k);
      C(iu, jv) = det(m);
    endfor
  endfor
endfunction