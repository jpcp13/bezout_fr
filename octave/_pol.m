function m = _pol(uu, vv, k)
  global n f df
  m = zeros(n, n);
  for j = 1:n
    if abs(uu(j) - vv(j)) < 1e-15
      for i = 1:n
        fuu = map(f{i}, [vv(1:j-1) uu(j:n)]);
        dfuu = map(df{i, j}, [vv(1:j-1) uu(j:n)]);
        if j != k
          m(i, j) = dfuu;
        else
          m(i, j) = uu(j)*dfuu - fuu;
        endif
      endfor
    else
      for i = 1:n
        fuu = map(f{i}, [vv(1:j-1) uu(j:n)]);
        fvv = map(f{i}, [vv(1:j) uu(j+1:n)]);
        if j != k
          m(i, j) = (fuu - fvv)/(uu(j) - vv(j));
        else
          m(i, j) = (vv(j)*fuu - uu(j)*fvv)/(uu(j) - vv(j));
        endif
      endfor
    endif
  endfor
endfunction