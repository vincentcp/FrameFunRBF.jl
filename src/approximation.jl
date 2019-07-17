
"""
Make a cartesian product equispaced grid on `[-1,1]^N`.
"""
function equispaced_grid(::Val{N}, n::Int) where {N}
    grid1 = PeriodicEquispacedGrid{T}(n)
    grid = ProductGrid( ntuple(t->grid1, Val{N}())... )
end

make_centers(grid::AbstractGrid) = collect(grid)

restrict_grid(grid, domain) = MaskedGrid(grid, domain)

function rbfapprox(rbf::RBF, centers, points, f; verbose = false, method = :svd)
    basis = RBFDictionary(centers, rbf)
    A = evaluation_matrix(basis, points)
    B = f.(points)
    if method == :LAPACK
        coef, rank = Base.LinAlg.LAPACK.gelsy!(A,B)
        if verbose
            println("Effective rank: ", rank)
        end
    elseif method == :svd
        threshold = 1e-14
        usv = LAPACK.gesdd!('S',A)
        s = usv[2]
        maxindex = findlast(s.>threshold)
        if verbose
            println("Effective svd-rank: ", maxindex)
        end
        sreg = 1 ./ s[1:maxindex]
        ureg = usv[1][:,1:maxindex]
        vreg = usv[3][1:maxindex,:]
        coef = vreg' * diagm(sreg) * (ureg'*B)
    else
        error("Please supply a valid method for rbfapprox")
    end

    Expansion(basis, coef)
end
