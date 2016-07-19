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

def _F(n, deg, fshape, m):
	F = [];
	for i in range(n):
		degrees = _degrees(n, deg, m)
		coeffs = np.random.randint(-10, 10, size=(m))
		f = sparse2prism(deg, fshape, degrees, coeffs)
		F.append(f)
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
