import matplotlib.pyplot as plt
import bezmat as bm
import scipy.io as sio
# import time
# import sys



deg = [2, 3, 4]
m = 20	# nombre de monomes
B = bm._bezout(deg, m)
B = bm.block_triang(deg, B)
plt.spy(B[0]); plt.grid(); plt.savefig('ref.png')

tB = bm._t(B)
bls = bm.block_size(deg)

sio.savemat('np_B.mat', {'B':tB, 'deg':deg, 'bls':bls, 'nb_monomials':m})
