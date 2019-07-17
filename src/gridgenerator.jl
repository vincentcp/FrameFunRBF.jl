
"""
A `GridGenerator` represents a family of grids.

The `generate` function maps a parameter value to a grid.
"""
abstract type GridGenerator{T} end

parameter(grid::GridGenerator) = 1


struct ModelGridGenerator{T} <: GridGenerator{T}
    modelgrid   ::  AbstractGrid{T}
end

parameter(gen::ModelGridGenerator) = dimensions(gen.modelgrid)

generate(gen::ModelGridGenerator, param) = resize(gen.modelgrid, param)


struct SubGridGenerator{T} <: GridGenerator{T}
    domain      ::  Domain{T}
    generator   ::  GridGenerator{T}
end

parameter(gen::SubGridGenerator) = parameter(gen.generator)

generate(gen::SubGridGenerator, param) = subgrid(generate(gen.generator, param), gen.domain)


struct EquispacedGrids{T} <: GridGenerator{T}
    domain  ::  Domain{T}
end

equispaced_grid(box::AbstractInterval{T}, n) where {T} =
    EquispacedGrid{T}(n, infimum(box), supremum(box))

equispaced_grid(box::ProductDomain, n) =
    ProductGrid(equispaced_grid.(elements(box), n)...)

generate(gen::EquispacedGrids, param) = equispaced_grid(gen.domain, param)


struct HexagonalGrids{T} <: GridGenerator{SVector{2,T}}
    bounding_box   ::  Vector{T}
end

function hexagonal_grid(N, bounding_box)
    Rad3Over2 = sqrt(3) / 2
    x = 0:1:round(sqrt(N)-1)
    X = [c for d in x, c in x]
    Y = [d for d in x, c in x]
    n = length(x)
    X = Rad3Over2 * X
    X /= maximum(X)
    Y[:,2:2:end] .+= 0.5
    Y /= maximum(Y)
    Lx = bounding_box[2]-bounding_box[1]
    Ly = bounding_box[4]-bounding_box[3]
    X = X[:] * Lx .+ bounding_box[1]
    Y = Y[:] * Ly .+ bounding_box[3]
    G = map(SVector, X, Y)
end

generate(gen::HexagonalGrids, param::Int) = hexagonal_grid(param, gen.bounding_box)
