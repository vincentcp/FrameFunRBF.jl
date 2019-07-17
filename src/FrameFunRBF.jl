module FrameFunRBF

using DomainSets, BasisFunctions, FrameFun
using StaticArrays, LinearAlgebra
using Printf
using CompactTranslatesDict

BF = BasisFunctions

import Base:
    length,
    size

import BasisFunctions:
    support,
    name,
    unsafe_eval_element,
    plotgrid

import FrameFun:
    dictionary,
    SamplingStyle,
    SolverStyle


export RBF, SmoothRBF, PiecewiseSmoothRBF
export Multiquadric,
    InverseMultiquadric,
    InverseQuadratic,
    Gaussian,
    PolyharmonicSpline,
    ThinPlateSpline

export RBFPlatform,
    RBFParLinear,
    RBFParConstant,
    RBFParSquareRoot,
    rbfselect

export RBFDictionary,
    PeriodicRBFDictionary


include("gridgenerator.jl")
include("rbf.jl")
include("rbf_dict.jl")
include("approximation.jl")
include("platform.jl")

end
