# module BezoutModule

C = Complex{Float64}
type mat_bool
    A::Array{Float64, 3}
    z::Bool
end
type mat_int
    A::Array{Float64, 3}
    r::Int
end

function giv(v::Array{Float64, 1})
    """ computes givens rotation """
    a, b = v[1], v[2]
    if abs(b) < epsi return eye(2) end
    if abs(b) > abs(a)
        t = -a/b; s = 1/sqrt(1 + t*t); c = s*t
    else
        t = -b/a; c = 1/sqrt(1 + t*t); s = c*t
    end
    return hcat([c, s], [-s, c])
end

function qrp(B::Array{Float64, 3}, idx::UnitRange{Int64} , idy::UnitRange{Int64} , k::Int,
    post_mult::Bool)
    """
    operates on left side of matrix - premultiplication by Q'
    """
    A = copy(B)
    a = A[idx, idy, k]
    F = qrfact(a, Val{true})
    Q, R, p = F[:Q], F[:R], F[:p]
    for l = 0:n
        A[idx, :, l+1] = Q'*A[idx, :, l+1]
        if post_mult
            A[:, idy, l+1] = A[:, idy, l+1][:, p]
        end
    end
    f = find(abs(diag(F[:R])).>epsi)
    r = isempty(f) ? 0 : maximum(f)
    return mat_int(A, r)
end

function ru(B::Array{Float64, 3}, i_max::Int, j_min::Int, j_max::Int)
    """
    operates on right side of matrix - postmultiplication by Tb
    """
    A = copy(B)
    ix = 1:i_max
    jt = j_min : j_min + i_max - 1
    jb = j_min + i_max : j_max
    T = A[ix, jt, n+1]
    b = A[ix, jb, n+1]
    Tb = triu(T)\b
    for l = 0:n
        F = A[:, jt, l+1]
        g = A[:, jb, l+1]
        A[:, j_min:j_max, l+1] = hcat(g - F*Tb, F)
    end
    return A
end

function triang_block(B::Array{Float64, 3})
    """ operates on both sides of matrix - calls qrp and ru """
    A = copy(B)
    loss = 0
    for k = 1:fn
        idx = (k-1)*s+1-loss : k*s
        idy = (k-1)*s+1 : k*s
        qrpA = qrp(A, idx, idy, n+1, true)
        A, r = qrpA.A, qrpA.r
        i_max = (k-1)*s+r-loss
        j_min = loss+1
        j_max = k*s
        A = ru(A, i_max, j_min, j_max)
        loss += s-r
    end
    return mat_int(A, loss)
end


function reflect(f::Function)
    function rf(B::Array{Float64, 3}, ij::Int, k::Int)
        ftB = f(permutedims(B, [2, 1, 3]), ij, k)
        if isa(ftB, mat_bool)
            return mat_bool(permutedims(ftB.A, [2, 1, 3]), ftB.z)
        elseif isa(ftB, mat_int)
            return mat_int(permutedims(ftB.A, [2, 1, 3]), ftB.r)
        else
            return permutedims(ftB, [2, 1, 3])
        end
        return rf
    end
end

function col_giv(B::Array{Float64, 3}, j::Int, k::Int)
    A = copy(B)
    Dx, Dy = size(A[:, :, n+1])
    for i = Dx-1:-1:1
        v = A[i:i+1, j, k]
        G = giv(v)
        for l = 0:n
            A[i:i+1, :, l+1] = G*A[i:i+1, :, l+1]
        end
    end
    return vcat(A[2:Dx, :, :], zeros(1, Dy, n+1))
end
row_giv(B::Array{Float64, 3}, j::Int, k::Int) = reflect(col_giv)(B, j, k)

function sweep_l(B::Array{Float64, 3}, m::Int, k::Int)
    A = copy(B)
    Dx, Dy = size(A)
    z = false
    for i = 1:m-1
        if abs(A[i, m-i+1, k]) < epsi
            v = A[i:i+1, m-i, k]
            G = giv(v)
            for l = 0:n
                A[i:i+1, :, l+1] = G*A[i:i+1, :, l+1]
            end
            z = true
        end
    end
    return mat_bool(A, z)
end
sweep_r(B::Array{Float64, 3}, m::Int, k::Int) = reflect(sweep_l)(B, m, k)
function sweep(A::Array{Float64, 3}, m, k)
    swAl = sweep_l(A, m, k)
    A, zl = swAl.A, swAl.z
    swAr = sweep_r(A, m, k)
    A, zr = swAr.A, swAr.z
    return mat_bool(A, zr | zl)
end

function h_relations(B::Array{Float64, 3}, m::Int, k::Int)
    A = copy(B)
    Dx, Dy = size(A[:, :, n+1])
    idx = m+1:Dx
    idy = 1:Dy
    return qrp(A, idx, idy, k, false)
end
v_relations(B::Array{Float64, 3}, i::Int, k::Int) = reflect(h_relations)(B, i, k)

function h_reduct(B::Array{Float64, 3}, m::Int, k::Int)
    hrelB = h_relations(B, m, k)
    B, r = hrelB.A, hrelB.r
    for i = m+1:m+r
        B = row_giv(B, i, k)
        swB = sweep(B, m, n+1)
        B, z = swB.A, swB.z
        if z
            m -= 1
        end
    end
    return mat_int(B, m)
end

v_reduct(B::Array{Float64, 3}, m::Int, k::Int) = reflect(h_reduct)(B, m, k)

function grand_reduct(B::Array{Float64, 3}, m::Int)
    for k = 1:n
        hredB = h_reduct(B, m, k)
        B, m = hredB.A, hredB.r
        vredB = v_reduct(B, m, k)
        B, m = vredB.A, vredB.r
    end
    return mat_int(B, m)
end

function plots(bez::Array{Float64, 3}, beztri::Array{Float64, 3}, B::Array{Float64, 3}, i_start::Int)
    clf()
    subplot(1, 2, 1)
    spy(bez[:, :, n+1])
    F = qrfact(bez[:,:,n+1], Val{true})
    Q, R, p = F[:Q], F[:R], F[:p]
    subplot(1, 2, 2)
    plot(log10(1e-16 + abs(diag(R))), "r*")
    savefig("bez_diagR.png")
    clf()
    subplot(1, 2, 1)
    spy(beztri[:, :, n+1])
    subplot(1, 2, 2)
    plot(log10(1e-16 + abs(diag(flipdim(B[:, :, n+1], 2), i_start))), "r*")
    savefig("beztri_diagR.png")
end
