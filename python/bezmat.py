#######################" fonctions python numpy scipy #########################
import numpy as np
import scipy.linalg as la
import scipy.fftpack as ft


def _GH(n, deg, dx, dy):
	fn = np.prod(range(1, n+1))
	Gx = []; Hx = []; Gy = []; Hy = []
	for j in range(n):
		dft = ft.fft(np.eye(2*fn*deg[j]))
		gx = dft[0::2*fn/(j+1), :deg[j]+1]
		gy = dft[1::2*fn/(n-j), :deg[j]+1]
		hx = dft[0::2*fn/(j+1), :dx[j]]
		hy = dft[1::2*fn/(n-j), :dy[j]]
		Gx.append(gx); Hx.append(hx); Gy.append(gy); Hy.append(hy)
	return Gx, Gy, Hx, Hy

def _HK(Hx, Hy):
	n = len(Hx)
	H = 1; K = 1
	for k in range(n):
		H = np.kron(Hx[k], H); K = np.kron(Hy[k], K)
	return H, K

def _degrees(n, deg, m):
	degrees = np.empty((m, n), dtype=int)
	for j in range(n):
		degrees[:, j] = np.random.randint(deg[j]+1, size=(m))
	return degrees

def sparse2prism(deg, fshape, degrees, coeffs):
	l = degrees.shape[0]
	f = np.zeros(np.asarray(deg)+1, dtype=int)
	for k in range(l):
		f[tuple(degrees[k])] = coeffs[k]
	return f

def _degs2latex(filename, f_degs, n, m, i):
	outfile = open(filename, 'a')
	outfile.write('$$\\begin{array}{c|' + 'c'*m + '} \n')
	outfile.write(('f_%1d ' + '&'*m + '\\\ \n') %(i+1))
	outfile.write('\\hline \n')
	for j in range(n):
		outfile.write(('d_%1d ') %(j+1))
		for k in range(m):
			outfile.write('&%2d' %f_degs[k, j])
		outfile.write('\\\ \n')
	outfile.write('\\end{array}$$ \n')
	outfile.close()

def _coeffs2latex(filename, f_coeffs, n, m):
	outfile = open(filename, 'w')
	outfile.write('$$\\begin{array}{' + 'r'*n + '} \n')
	for i in range(n-1):
		outfile.write('f_%1d & ' %(i+1))
	outfile.write('f_%1d \\\ \n' %n)
	outfile.write('\\hline \n')
	for k in range(m):
		for i in range(n-1):
			outfile.write('%2d & ' %f_coeffs[k, i])
		outfile.write('%2d \\\ \n' %f_coeffs[k, n-1])
	outfile.write('\\end{array}$$ \n')
	outfile.close()

def _F(n, deg, fshape, m):
	F = [];
	outfile = open('../tex/degrees.tex', 'w')
	outfile.close()
	coeffs_table = np.empty((m, n), dtype=int)
	for i in range(n):
		degrees = _degrees(n, deg, m)
		_degs2latex('../tex/degrees.tex', degrees, n, m, i)
		coeffs = np.random.randint(-2, 2, size=(m))
		coeffs_table[:, i] = coeffs
		f = sparse2prism(deg, fshape, degrees, coeffs)
		F.append(f)
	_coeffs2latex('../tex/coeffs.tex', coeffs_table, n, m)
	degrees = np.zeros((1, n), dtype=int)
	coeffs = [1]
	f = sparse2prism(deg, fshape, degrees, coeffs)
	F.append(f)
	for i in range(n):
		degrees = np.zeros((1, n), dtype=int)
		degrees[0, i] = 1
		coeffs = [1]
		f = sparse2prism(deg, fshape, degrees, coeffs)
		F.append(f)
	return F

def evaluate(F, n, fshape, dx, dy, Gx, Gy, i, j):
	f = F[i].transpose(range(n-1, -1, -1)).reshape(np.prod(fshape[j:],dtype=int), np.prod(fshape[:j], dtype=int))
	L = 1; R = 1
	for k in range(j):
		o = np.ones((dx[k], 1))
		L = np.kron(o, L); R = np.kron(Gy[k], R)
	for k in range(j, n):
		o = np.ones((dy[k], 1))
		R = np.kron(o, R); L = np.kron(Gx[k], L)
	return L.dot(f).dot(R.T)

def _J(F, n, fshape, dx, dy, Gx, Gy):
	Dx, Dy = np.prod(dx), np.prod(dy)
	J = np.empty((Dx, Dy, 2*n+1, n+1), dtype=complex)
	for i in range(2*n+1):
		for j in range(n+1):
			J[:, :, i, j] = evaluate(F, n, fshape, dx, dy, Gx, Gy, i, j)
	return J

def _C(n, J):
	Dx, Dy = J.shape[:2]
	C = np.empty((n+1, Dx, Dy), dtype=complex)
	indices_vol = range(n, 2*n+1)
	for i in range(n+1):
		indices = range(n) + [n+i]
		vol = np.linalg.det(J[:, :, indices_vol, :])
		num = np.linalg.det(J[:, :, indices, :])
		C[i] = num/vol
	return C

def _B(n, C, H, K):
	Dx, Dy = C.shape[-2:]
	B = np.empty((n+1, Dx, Dy), dtype=complex)
	for i in range(n+1):
		HC  = np.conjugate(H.T).dot(C[i])/Dx
		KHC = np.conjugate(K.T).dot(HC.T)/Dy
		B[i] = KHC.T
	return np.around(B).real.astype(int)

def _bezout(deg, m):
	fshape = [d+1 for d in deg]
	n = len(deg)
	dx = [(i+1)*deg[i] for i in range(n)]
	dy = [(n-i)*deg[i] for i in range(n)]
	Gx, Gy, Hx, Hy = _GH(n, deg, dx, dy)
	H, K = _HK(Hx, Hy)
	F = _F(n, deg, fshape, m)
	J = _J(F, n, fshape, dx, dy, Gx, Gy)
	C = _C(n, J)
	B = _B(n, C, H, K)
	return B

def build_ixy(deg, k):
	n = len(deg)
	ix = ()
	iy = ()
	for j in range(n):
		if j < k:
			ix = ix + (slice(None),)
			iy = (slice(None),) + iy
		elif j == k:
			ix = ix + (slice(-deg[k], None),)
			iy = (slice(-deg[n-1-k], None),) + iy
		else:
			ix = ix + (slice(None, -deg[j]),)
			iy = (slice(None, -deg[n-1-j]),) + iy
	return ix, iy

def permut(deg):
	n = len(deg)
	dx = [(i+1)*deg[i] for i in range(n)]
	dy = [(n-i)*deg[i] for i in range(n)]
	Dx, Dy = np.prod(dx), np.prod(dy)
	aax = np.arange(Dx, dtype=int).reshape(dx[::-1]).transpose(range(n-1,-1,-1))
	aay = np.arange(Dy, dtype=int).reshape(dy[::-1]).transpose(range(n-1,-1,-1))
	ax = np.zeros(0, dtype=int)
	ay = np.zeros(0, dtype=int)
	for k in range(n):
		ix, iy = build_ixy(deg, k)
		tx = (np.arange(n)*(n-1) + k) % n
		ty = (n-1 - tx) % n
		ay = np.concatenate(( ay, np.transpose(aay[iy], ty).reshape(Dy/n) ))
		ax = np.concatenate(( ax, np.transpose(aax[ix], tx).reshape(Dx/n) ))
	return ax, ay

def block_triang(deg, B):
	Btri = np.copy(B)
	n = len(deg)
	ax, ay = permut(deg)
	for k in range(n+1):
		Btri[k] = np.fliplr(B[k][np.ix_(ax,ay)])
	return Btri

def block_size(deg):
	n = len(deg)
	dx = [(i+1)*deg[i] for i in range(n)]
	Dx = np.prod(dx)
	bls = np.zeros(0, dtype=int)
	for i in range(n):
		bls = np.concatenate(( bls, np.ones(deg[i], dtype=int)*Dx/(n*deg[i]) ))
	return bls

def _save2text(B):
	x = B.reshape(np.prod(B.shape))
	np.savetxt('../julia/B.txt', x, fmt='%8d')
