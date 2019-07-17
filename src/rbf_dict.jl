
struct RBFDictionary{S,T,R <: RBF{T}} <: Dictionary{S,T}
    centers ::  Vector{S}
    rbf     ::  R
    domain  ::  Domain{S}
end

RBFDictionary(centers::Vector{S}, rbf::RBF) where {S} = RBFDictionary(centers, rbf, DomainSets.FullSpace{S}())

size(dict::RBFDictionary) = size(dict.centers)

domain(dict::RBFDictionary) = dict.domain

name(dict::RBFDictionary) = "A list of radial basis functions ($(typeof(dict.rbf)))"

center(dict::RBFDictionary, idx) = dict.centers[idx]

support(dict::RBFDictionary) = dict.domain

unsafe_eval_element(dict::RBFDictionary, idx, x) =
    dict.rbf(norm(x - center(dict,idx)))

plotgrid(dict::RBFDictionary, n) = equispaced_grid(boundingbox(dict.domain), n)


import CompactTranslatesDict:
    CompactTranslationDict,
    eval_kernel,
    kernel_span

struct PeriodicRBFDictionary{T,R <: RBF{T}} <: CompactTranslatesDict.CompactTranslationDict{T}
    rbf     ::  R
    n       ::  Int
end

eval_kernel(dict::PeriodicRBFDictionary, x) = dict.rbf(x)

function kernel_span(dict::PeriodicRBFDictionary)
    D = approximate_support(dict.rbf)
    -D..D
end
