import matplotlib.pyplot as plt
import bezmat as bm
import scipy.io as sio
# import time
# import sys

deg = [2,2,2,2]
m = 20	# nombre de monomes
B = bm._bezout(deg, m)
B = bm.block_triang(deg, B)

bm._save2text(B)

plt.spy(B[0]); plt.grid(); plt.savefig('../png/ref.png')
