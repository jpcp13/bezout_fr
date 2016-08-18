#######################" fonctions python numpy scipy #########################

def _GH():
	Gx = []; Hx = []; Gy = []; Hy = []
	for j in range(n):
		dft = ft.fft(np.eye(2*fn*deg[j]))
		gx = dft[0::2*fn/(j+1), :deg[j]+1]
		gy = dft[1::2*fn/(n-j), :deg[j]+1]
		hx = dft[0::2*fn/(j+1), :dx[j]]
		hy = dft[1::2*fn/(n-j), :dy[j]]
		Gx.append(gx); Hx.append(hx); Gy.append(gy); Hy.append(hy)
	return Gx, Gy, Hx, Hy

def _HK():
	H = 1; K = 1
	for k in range(n):
		H = np.kron(Hx[k], H); K = np.kron(Hy[k], K)
	return H, K

def evaluate(i, j):
	f = F[i].transpose(range(n-1, -1, -1)).reshape(prod(fshape[j:]), prod(fshape[:j]))
	L = 1; R = 1
	for k in range(j):
		o = np.ones((dx[k], 1))
		L = np.kron(o, L); R = np.kron(Gy[k], R)
	for k in range(j, n):
		o = np.ones((dy[k], 1))
		R = np.kron(o, R); L = np.kron(Gx[k], L)
	return L.dot(f).dot(R.T)

def _J():
	J = np.empty((Dx, Dy, 2*n+1, n+1), dtype=complex)
	for i in range(2*n+1):
		for j in range(n+1):
			J[:, :, i, j] = evaluate(i, j)
	return J

def _C():
	C = np.empty((n+1, Dx, Dy), dtype=complex)
	indices_vol = range(n, 2*n+1)
	for i in range(n+1):
		indices = range(n) + [n+i]
		vol = np.linalg.det(J[:, :, indices_vol, :])
		num = np.linalg.det(J[:, :, indices, :])
		C[i] = num/vol
	return C

def _B():
	B = np.empty((n+1, Dx, Dy), dtype=complex)
	for i in range(n+1):
		HC  = conjugate(H.T).dot(C[i])/Dx
		KHC = conjugate(K.T).dot(HC.T)/Dy
		B[i] = KHC.T
	return np.around(B).real.astype(int)

def build_ixy(k):
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

def permut():
	aax = np.arange(Dx, dtype=int).reshape(dx[::-1]).transpose(range(n-1,-1,-1))
	aay = np.arange(Dy, dtype=int).reshape(dy[::-1]).transpose(range(n-1,-1,-1))
	ax = np.zeros(0, dtype=int)
	ay = np.zeros(0, dtype=int)
	for k in range(n):
		ix, iy = build_ixy(k)
		tx = (np.arange(n)*(n-1) + k) % n
		ty = (n-1 - tx) % n
		ay = np.concatenate(( ay, np.transpose(aay[iy], ty).reshape(Dy/n) ))
		ax = np.concatenate(( ax, np.transpose(aax[ix], tx).reshape(Dx/n) ))
	return ax, ay

def block_triang():
	ax, ay = permut()
	for k in range(n+1):
		B[k] = np.fliplr(B[k][np.ix_(ax,ay)])
	return B

def block_size():
	bls = np.zeros(0, dtype=int)
	for i in range(n):
		bls = np.concatenate(( bls, np.ones(deg[i], dtype=int)*Dx/(n*deg[i]) ))
	return bls

############################################# fcts sage #######################################

def rand_poly(j, m):
	# construit un polynome random en j+1 variables 0..j, de degres max(deg, m)
	p = 0
	for k in range(min(deg[j], m) + 1):
		if j > 0:
			p += rand_poly(j-1, m-k)*x[j]**k
		else:
			coeff = ZZ.random_element(0.5, distribution='gaussian')
			p += coeff*x[j]**k
	return p

def poly2prism(p):
	z = np.zeros(fshape)
	mns = p.monomials(); cfs = p.coefficients()
	for k in range(len(mns)):
		z[mns[k].degrees()] = cfs[k]
	return z

def poly2sparse(P):
	degrees = np.zeros((0, n))
	coeffs = np.zeros((0))
	nb_monomials = np.zeros(n, dtype=int)
	for j in range(n):
		p = P[j]
		mns = p.monomials();
		degs = np.zeros((len(mns), n))
		for k in range(len(mns)):
			degs[k, :] = np.asarray(mns[k].degrees())
		cfs = np.asarray(p.coefficients())
		degrees = np.concatenate((degrees, degs))
		nb_monomials[j] = len(mns)
		coeffs = np.concatenate((coeffs, cfs))
	return degrees, coeffs, nb_monomials

def sparse2prism(sp):
	degrees = sp[0]
	l = degrees.shape[0]
	coeffs = sp[1]
	f = np.zeros(np.asarray(deg)+1, dtype=float)
	for k in range(l):
		f[tuple(degrees[k])] = coeffs[k]
	return f

def shrink(axe):
	iz = np.all(B == 0, axis=axe)
	pix = np.prod(iz, axis=0)
	ix = np.nonzero(1 - pix)[0]
	return list(ix)

def block_right_kernel(B):
	Dx, Dy = B.nrows(), B.ncols()
	s = bls[0]
	X1 = matrix(QQ, 0, 0)
	for k in range(sum(deg)):
		B0 = B[Dx-k*s-s:Dx-k*s, Dy-k*s-s:Dy-k*s]
		B1 = B[Dx-k*s-s:Dx-k*s, Dy-k*s:Dy]
		r1 = X1.ncols()
		b1x1 = B1*X1
		B0a = B0.augment(b1x1)
		Y0 = B0a.right_kernel().basis_matrix().transpose()
		X0 = Y0[0:s, :]
		X1 = X1*Y0[s:s+r1:, :]
	 	X1 = X0.transpose().augment(X1.transpose()).transpose()
	return X1

def reduct():
	nr, nc = BB[0].nrows(), BB[0].ncols()
	ker_basis = BB[0].left_kernel().basis_matrix()
	print 'ker_size = ', ker_basis.nrows()
	relations = matrix(QQ, 0 , nc)
	for k in range(1, n+1):
		new_relat = ker_basis*BB[k]
		relations = block_matrix(2, 1, [relations, new_relat])
	r = rank(relations)
	print 'r = ', r
	rkbt = relations.right_kernel().basis_matrix().transpose()
	BB_out = []
	for k in range(n+1):
		BB_out.append(BB[k]*rkbt)
	den = BB_out[0].denominator()
	print 'den = ', den
	return r, BB_out

def _Bred():
	nr, nc = BB[0].nrows(), BB[0].ncols()
	d = rank(BB[0])
	print 'd =', d
	piv_cols, piv_rows = BB[0].pivots(), BB[0].transpose().pivots()
	Bred = np.empty((n+1, d, d), dtype=float)
	for j in range(n+1):
		Bred[j] = BB[j][piv_rows, piv_cols] + matrix(RR, d, d)
	return Bred

def _XX_chow():
	nr, nc = BB[0].nrows(), BB[0].ncols()
	d = rank(BB[0])
	print 'd =', d
	piv_cols, piv_rows = BB[0].pivots(), BB[0].transpose().pivots()
	nBB = np.empty((n+1, d, d), dtype=float)
	XX = np.empty((n, d, d), dtype=float)
	for k in range(n+1):
		nBB[k, :, :] = np.array(BB[k][piv_rows, piv_cols])
	for k in range(n):
		XX[k, :, :] = la.solve(nBB[0], nBB[k+1])
	chow_mat = np.zeros((d, d))
	for k in range(n):
		chow_mat += np.random.random()*nBB[k+1]
	Ex = la.eig(chow_mat, nBB[0])
	Wx = Ex[1]
	X = np.empty((n, d), dtype=complex)
	for k in range(n):
		X[k] = la.solve(Wx, XX[k].dot(Wx)).diagonal(offset=0)
	return XX, X


def _jPZ(X):
	d = X.shape[1]
	jz = np.empty((n, n), dtype=complex)
	Pz = np.empty((n, 1), dtype=complex)
	jPZ = np.empty((n, d), dtype=complex)
	for k in range(d):
		z = tuple(X[:, k])
		for i in range(n):
			Pz[i] = P[i](z)
#			for j in range(n):
#				jz[i, j] = jac[i, j](z)
#		if abs(la.det(jz)) > 1e-4:
#			jPZ[:, k:k+1] = la.solve(jz, Pz)
#		else:
		jPZ[:, k:k+1] = Pz
	return jPZ

def plot_octave():
	os.system("sh test.sh")
	rac_real = np.loadtxt('../txt/rac_real.txt')
	rac_imag = np.loadtxt('../txt/rac_imag.txt')
	rac = rac_real + rac_imag*1j
	X1 = rac.transpose()
	jPZ = _jPZ(rac.transpose())
	plt.subplot(1, 2, 2)
	#plt.plot(np.log10(2**-52 + abs(np.transpose(jPZ[:,:]))))
	for i in range(n):
		hh = plt.hist(np.log10(2**-52 + abs(jPZ[i])), 50)
	plt.grid();plt.savefig('../png/ref.png')

def int2tex(filename, i):
	outfile = open(filename, 'w')
	outfile.write('$%3d$' %i)
	outfile.close()

def list2tex(filename, l):
	outfile = open(filename, 'w')
	outfile.write('[')
	for i in range(len(l)-1):
		li = l[i]
		outfile.write('%1d ,' %li)
	l1 = l[-1]
	outfile.write('%1d]' %l1)
	outfile.close()

def _rank(B0):
	nr, nc = B0.nrows(), B0.ncols()
	ker_basis = B0.left_kernel().basis_matrix()
	return Dx - ker_basis.nrows()
##################" debut programme python / sage ######################

import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as la
import scipy.fftpack as ft
import time
import sys
import scipy.io as sio

deg = [2,2,2]
list2tex('../txt/deg.txt', deg)

m = 15
n = len(deg)
int2tex('../txt/n.txt', n)

R = PolynomialRing(QQ, 'x', n)
x = R.gens()
xx = [x[0]**0] + list(x)

fshape = [d+1 for d in deg]
dx = [(i+1)*deg[i] for i in range(n)]
dy = [(n-i)*deg[i] for i in range(n)]
fn, Dx, Dy = factorial(n), prod(dx), prod(dy)
int2tex('../txt/Dx.txt', Dx)

P = [rand_poly(n-1, m) for i in range(n)] + xx
#~ a = [ZZ.random_element() for _ in range(7)]
#~ P = [a[0]*x[0]**3*x[1]**2 + a[1]*x[0] + a[2]*x[1]**2 + a[3], a[4]*x[0]*x[1]**4 + a[5]*x[0]**3 + a[6]*x[1]] + xx

# P = P + xx
# degrees, coeffs, nb_monomials = poly2sparse(P)

jac = jacobian(P[:n], x)
F = [poly2prism(p) for p in P]

Gx, Gy, Hx, Hy = _GH()
H, K = _HK()
J = _J()
C = _C()

B = _B()
plt.close()
#plt.subplot(1, 2, 1)
plt.spy(B[0])
plt.savefig('../png/bez.png')
qq, rr, pp = la.qr(B[0], pivoting=True)
dr = rr.diagonal()
plt.close()
plt.semilogy(abs(dr)+1e-18, '*')
plt.grid()
plt.savefig('../png/bez_diag.png')

bls = block_size()

B = block_triang()
plt.close()
plt.spy(B[0])
plt.savefig('../png/beztri.png')
DR = np.empty(Dx)
idx = 0
for s in bls:
	bb = B[0][idx:idx+s, idx:idx+s]
	qq, rr, pp = la.qr(bb, pivoting=True)
	dr = rr.diagonal()
	DR[idx:idx+s] = dr
	idx += s
plt.close()
plt.semilogy(abs(DR)+1e-18, '*')
plt.grid()
plt.savefig('../png/beztri_diag.png')





print 'debut sage'
BB = []
for k in range(n+1):
	Bk = matrix(QQ, B[k])
	BB.append(Bk[:, :])

int2tex('../txt/dim0.txt', _rank(BB[0]))

r = 1
while r > 0:
	r, BB = reduct()

Bred = _Bred()



sio.savemat('np_B.mat', {'Bred':np.transpose(Bred.astype(float), (1, 2, 0)), 'B':np.transpose(B.astype(float), (1, 2, 0)), 'deg':deg, 'bls':bls,})

# 'degrees':degrees,
# 'coeffs':coeffs,
# 'nb_monomials':nb_monomials

XX, X = _XX_chow()
jPZ = _jPZ(X)

plt.close()
plt.subplot(1, 2, 1)
#plt.plot(np.log10(2**-52 + abs(np.transpose(jPZ[:,:]))))
for i in range(n):
	hh = plt.hist(np.log10(2**-52 + abs(jPZ[i])), 50)
plt.grid();plt.savefig('../png/ref.png')

plot_octave()

I = R.ideal(P[:n])
dim = I.vector_space_dimension()
print 'dim =', dim

int2tex('../txt/dim.txt', dim)
