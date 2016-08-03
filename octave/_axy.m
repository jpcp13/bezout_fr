function [ax2, ay2] = _axy()
  global n D degx degy
  ax = reshape(1:D, fliplr(degx));
  ax = permute(ax, n:-1:1);
  ax = reshape(ax, 1, D);
  ay = reshape(1:D, fliplr(degy));
  ay = permute(ay, n:-1:1);
  ay = reshape(ay, 1, D);
  aax = reshape(ax, degx);
  aay = reshape(ay, degy);
  ax2 = []; ay2 = [];
  for k = 1:n
    [ix, iy] = _ixy(k);
    tx = mod( (0:n-1)*(n-1) + k-1, n) + 1;
    ty = mod(n-tx, n) + 1;
    aix = aax(ix{:});
    pix = permute(aix, fliplr(tx));
    ax2 = [ax2 reshape(pix, 1, D/n)];
    aiy = aay(iy{:});
    piy = permute(aiy, fliplr(ty));
    ay2 = [ay2 reshape(piy, 1, D/n)];
  endfor
endfunction
