import matplotlib.pyplot as plt
import bezmat as bm
import scipy.io as sio
# import time
# import sys

deg = [2,2,2,2]
m = 15	# nombre de monomes

B = bm._bezout(deg, m)
bm._save2text(B, '../julia/B.txt')

Btri = bm.block_triang(deg, B)
bm._save2text(Btri, '../julia/Btri.txt')

plt.subplot(1, 2, 1)
plt.spy(B[0]); plt.grid()
plt.subplot(1, 2, 2)
plt.spy(B[0]); plt.grid()
plt.savefig('../png/B0_B0tri.png')


