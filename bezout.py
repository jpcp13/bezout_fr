

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

def shrink(axe):
	iz = np.all(B == 0, axis=axe)
	pix = np.prod(iz, axis=0)
	ix = np.nonzero(1 - pix)[0]
	return list(ix)

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


def _jPZ():
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

##################" debut programme python / sage ######################

# import matplotlib.pyplot as plt
# import numpy as np
# import scipy.linalg as la
# import scipy.fftpack as ft
# import time
# import sys
# import scipy.io as sio
import bezmat as bm

deg = [2, 3]
fshape = [d+1 for d in deg]
n = len(deg)
dx = [(i+1)*deg[i] for i in range(n)]
dy = [(n-i)*deg[i] for i in range(n)]

m = 5	# nombre de monomes
# R = PolynomialRing(QQ, 'x', n)
# x = R.gens()
# xx = [x[0]**0] + list(x)


# P = [rand_poly(n-1, m) for i in range(n)] + xx
#~ a = [ZZ.random_element() for _ in range(7)]
#~ P = [a[0]*x[0]**3*x[1]**2 + a[1]*x[0] + a[2]*x[1]**2 + a[3], a[4]*x[0]*x[1]**4 + a[5]*x[0]**3 + a[6]*x[1]] + xx
# degrees, coeffs, nb_monomials = poly2sparse(P)

# jac = jacobian(P[:n], x)
# F = [poly2prism(p) for p in P]

Gx, Gy, Hx, Hy = bm._GH(n, deg, dx, dy)
H, K = bm._HK(Hx, Hy)
F = bm._F(n, deg, fshape, m)
J = bm._J(F, n, fshape, dx, dy, Gx, Gy)
C = bm._C(n, J)
B = bm._B(n, C, H, K)
#
#
# plt.spy(B[0]); plt.grid(); plt.savefig('ref.png')
#
# bls = block_size()
# B = block_triang()
#
#
# print 'debut sage'
# BB = []
# for k in range(n+1):
# 	Bk = matrix(QQ, B[k])
# 	BB.append(Bk[:, :])
# r = 1
# while r > 0:
# 	r, BB = reduct()
#
# #Bred = _Bred()
# #'Bred':np.transpose(Bred.astype(float), (1, 2, 0)),
#
#
#
# sio.savemat('np_B.mat', {'B':np.transpose(B.astype(float), (1, 2, 0)), 'deg':deg, 'bls':bls, 'degrees':degrees, 'coeffs':coeffs, 'nb_monomials':nb_monomials})
#
# #XX, X = _XX_chow()
#
# #jPZ = _jPZ()
#
# #plt.close()
# ##plt.plot(np.log10(2**-52 + abs(np.transpose(jPZ[:,:]))))
# #for i in range(n):
# #	hh = plt.hist(np.log10(2**-52 + abs(jPZ[i])), 50)
# #plt.grid();plt.savefig('ref.png')
#
# I = R.ideal(P[:n])
# dim = I.vector_space_dimension()
# print 'dim =', dim
#
