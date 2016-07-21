import numpy as np

m = 10
n = 3

coeffs = np.random.randint(-2, 2, size=(m, n))

outfile = open('coeffs.txt', 'w')

outfile.write('$$\\begin{array}{' + 'r'*n + '} \n')

for j in range(n-1):
	outfile.write('f_%1d & ' %(j+1))
outfile.write('f_%1d \\\ \n' %n)
outfile.write('\\hline \n')

for k in range(m):
	for j in range(n-1):
		outfile.write('%2d & ' %coeffs[k, j])
	outfile.write('%2d \\\ \n' %coeffs[k, n-1])

outfile.write('\\end{array}$$ \n')
outfile.close()
