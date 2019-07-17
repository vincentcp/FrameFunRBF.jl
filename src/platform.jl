
abstract type RBFParameterSelection end

struct RBFParLinear <: RBFParameterSelection
    c
end
struct RBFParSquareRoot <: RBFParameterSelection
    c
end
struct RBFParConstant <: RBFParameterSelection
    c
end

RBFparameter(p::RBFParLinear, n::Int) = p.c*n
RBFparameter(p::RBFParSquareRoot, n::Int) = p.c*sqrt(convert(typeof(p.c),n))
RBFparameter(p::RBFParConstant, n::Int) = p.c

name(p::RBFParLinear) = "$(@sprintf("%.2f", p.c)) * N"
name(p::RBFParSquareRoot) = "$(@sprintf("%.2f", p.c)) * sqrt(N)"
name(p::RBFParConstant) = "$(@sprintf("%.2f", p.c))"

rbfselect(rbf::RBF, n, select) = similar(rbf, RBFparameter(select, n))


struct RBFPlatform{R <: RBF} <: FrameFun.Platform
    rbf     ::  R
    domain  ::  Domain
    box     ::  Domain
    parameterstyle
    gridgen
end

default_selection(rbf::SmoothRBF{T}) where {T} = RBFParLinear(one(T))
default_selection(rbf::PiecewiseSmoothRBF) = RBFParConstant(rbf.p)

default_gridgenerator(rbf::RBF, domain::Domain{T}, box::Domain{T}) where {T} =
    EquispacedGrids{T}(box)

RBFPlatform(rbf::SmoothRBF{T}; domain = zero(T)..one(T),
        box = 2*boundingbox(domain),
        select = default_select(rbf),
        gridgen = default_gridgenerator(rbf, domain, box)) where {T} =
    RBFPlatform(rbf, domain, box, select, gridgen)

dictionary(pf::RBFPlatform, ppar) = dictionary(pf, ppar, pf.rbf, pf.domain, pf.box, pf.parameterstyle, pf.gridgen)

function dictionary(pf::RBFPlatform, n, rbf, domain, box, paramstyle, gridgen)
    rbfargs = RBFparameter(paramstyle, n)
    # centers = equispaced_centers(pf.domain, pf.box, n)
    centers = collect(generate(gridgen, n))
    RBFDictionary(centers[:], similar(rbf, rbfargs...), domain)
end

SamplingStyle(pf::RBFPlatform) = OversamplingStyle()
SolverStyle(pf::RBFPlatform) = DirectStyle()

# function equispaced_centers(domain::Domain, box, n)
#     g = equispaced_grid(box, n)
#     collect(g)
# end
#
# equispaced_grid(box::AbstractInterval{T}, n) where {T} =
#     EquispacedGrid{T}(n, infimum(box), supremum(box))
#
# equispaced_grid(box::ProductDomain, n) =
#     ProductGrid(equispaced_grid.(elements(box), n)...)
#
# periodic_equispaced_grid(box::AbstractInterval{T}, n) where {T} =
#     EquispacedGrid{T}(n, infimum(box), supremum(box))
#
# periodic_equispaced_grid(box::ProductDomain, n) =
#     ProductGrid(map(periodic_equispaced_grid, elements(box), n)...)
