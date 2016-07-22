
using PyPlot
include("bezoutFunctions.jl")

""" début programme """

epsi = 1e-9
n = 4
d = 2
D = prod(d:d:n*d)
fn = n*d # nombre de blocs
s = div(D, fn) # taille des blocs

b = readdlm("B.txt")
btri = readdlm("Btri.txt")
bb = reshape(b, D, D, n+1)
bbtri = reshape(btri, D, D, n+1)
bez = similar(bb)
beztri = similar(bbtri)
for k = 1:n
    bez[:, :, k] = transpose(bb[:, :, k+1])
    beztri[:, :, k] = transpose(bbtri[:, :, k+1])
    # bez[:, :, k] = flipdim(flipdim(bb[:, :, k+1], 1), 2)
    # beztri[:, :, k] = flipdim(flipdim(bbtri[:, :, k+1], 1), 2)
end
bez[:, :, n+1] = transpose(bb[:, :, 1])
beztri[:, :, n+1] = transpose(bbtri[:, :, 1])
# bez[:, :, n+1] = flipdim(flipdim(bb[:, :, 1], 1), 2)
# beztri[:, :, n+1] = flipdim(flipdim(bbtri[:, :, 1], 1), 2)


B = beztri[:, :, :]
tribloB = triang_block(B)
B, i_start = tribloB.A, tribloB.r
B = flipdim(B, 2)
# spy(abs(B[:, :, n+1]).>epsi)

m = size(B, 1) - i_start
println(m)


grB = grand_reduct(B, m)
B, m = grB.A, grB.r
println(m)

i_start = size(B, 1)-m
plot(log10(1e-16 + abs(diag(flipdim(B[:, :, n+1], 2), i_start))), "r*")
