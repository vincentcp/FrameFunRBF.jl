
"Any radial basis function."
abstract type RBF{T}
end


"A smooth radial basis function with a shape parameter `ε`."
abstract type SmoothRBF{T} <: RBF{T}
end

shape_parameter(rbf::SmoothRBF) = rbf.ε


"A piecewise smooth radial basis function."
abstract type PiecewiseSmoothRBF{T} <: RBF{T}
end

degree(rbf::PiecewiseSmoothRBF) = rbf.p


struct Multiquadric{T} <: SmoothRBF{T}
    ε   :: T
end

const MQ{T} = Multiquadric{T}

(rbf::Multiquadric)(r) = sqrt(1 + (rbf.ε*r)^2)

eval_derivative1(rbf::Multiquadric, r) = 1/rbf(r)*rbf.ε*r

similar(r::Multiquadric{T}, ε) where {T} = Multiquadric{T}(ε)


struct InverseMultiquadric{T} <: SmoothRBF{T}
    ε   :: T
end

const IMQ{T} = InverseMultiquadric{T}

(rbf::InverseMultiquadric)(r) = 1/sqrt(1 + (rbf.ε*r)^2)

eval_derivative1(rbf::InverseMultiquadric{T}, r) where {T} = (-one(T)/2)*rbf(r)^3*2*rbf.ε*r

similar(r::InverseMultiquadric{T}, ε) where {T} = InverseMultiquadric{T}(ε)


struct InverseQuadratic{T} <: SmoothRBF{T}
    ε   :: T
end

const IQ{T} = InverseQuadratic{T}

(rbf::InverseQuadratic)(r) = 1/(1 + (rbf.ε*r)^2)

similar(r::InverseQuadratic{T}, ε) where {T} = InverseQuadratic{T}(ε)


struct Gaussian{T} <: SmoothRBF{T}
    ε   :: T
end

const GA{T} = Gaussian{T}

(rbf::Gaussian)(r) = exp(-(rbf.ε*r)^2)

eval_derivative1(rbf::Gaussian, r) = -2*rbf.ε*r * rbf(r)

eval_derivative2(rbf::Gaussian, r) = -2*rbf.ε * rbf(r) - 2*rbf.ε*r*eval_derivative1(rbf,r)

similar(r::Gaussian{T}, ε) where {T} = Gaussian{T}(ε)

approximate_support(rbf::Gaussian{T}, threshold = eps(T)) where {T} = sqrt(-log(threshold))/rbf.ε


struct PolyharmonicSpline{T} <: PiecewiseSmoothRBF{T}
    p   ::  Int
end

const SPH{T} = PolyharmonicSpline{T}

PolyharmonicSpline(p::Int) = PolyharmonicSpline{Float64}(p)

(rbf::PolyharmonicSpline)(r) = r^(2*rbf.p+1)

similar(r::PolyharmonicSpline{T}, p) where {T} = PolyharmonicSpline{T}(p)


struct ThinPlateSpline{T} <: PiecewiseSmoothRBF{T}
    p   ::  Int
end

const TPS{T} = ThinPlateSpline{T}

ThinPlateSpline(p::Int) = ThinPlateSpline{Float64}(p)

(rbf::ThinPlateSpline)(r) = r^(2*rbf.p) * log(r)

similar(r::ThinPlateSpline{T}, p) where {T} = ThinPlateSpline{T}(p)
