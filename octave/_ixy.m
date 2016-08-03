function [ix, iy] = _ixy(k)
  global n deg degx degy
  ix = cell(1, n); iy = cell(1, n);
  for j = 1: k-1
    ix{j} = 1: degx(j);
    iy{n-j+1} = 1: degy(n-j+1);
  end
  ix{k} = degx(k)-deg(k)+1: degx(k);
  iy{n-k+1} = degy(n-k+1)-deg(n-k+1)+1: degy(n-k+1);
  for j = k+1: n
    ix{j} = 1: degx(j)-deg(j);
    iy{n-j+1} = 1: degy(n-j+1)-deg(n-j+1);
  end
endfunction
