import numpy as np

m = 10
n = 3
deg = [2, 2, 2]
degrees = np.empty((m, n), dtype=int)
for j in range(n):
	degrees[:, j] = np.random.randint(deg[j]+1, size=(m))

i = 0

outfile = open('degrees.txt', 'a')
outfile.write('$$\\begin{array}{c|' + 'c'*m + '} \n')
outfile.write(('f_%1d ' + '& '*m + '\\\ \n') %(i+1))
outfile.write('\\hline \n')
for j in range(n):
	outfile.write(('d_%1d ') %(j+1))
	for k in range(m):
		outfile.write('& %1d ' %degrees[k, j])
	outfile.write('\\\ \n')
outfile.write('\\end{array}$$ \n')
outfile.close()
