module MeixnerPollaczekPolynomials

import Base: axes, getindex, oneto
import ClassicalOrthogonalPolynomials: grammatrix, Inclusion, jacobimatrix, orthogonalityweight, OrthogonalPolynomial, Weight, ℝ, ∞
import LazyBandedMatrices: Tridiagonal
import RationalVals: RationalVal
import SpecialFunctions: gamma

export MeixnerPollaczek

struct MeixnerPollaczekWeight{T,L,P} <: Weight{T}
    λ::L
    ϕ::P
end
MeixnerPollaczekWeight(λ::L, ϕ::P) where {L,P} = MeixnerPollaczekWeight{promote_type(L, P)}(λ, ϕ)

axes(::MeixnerPollaczekWeight{T}) where {T} = (Inclusion{T}(ℝ),)
function getindex(w::MeixnerPollaczekWeight, x::Number)
    x ∈ axes(w, 1) || throw(BoundsError())
    abs2(gamma(w.λ + x * im)) * exp(x * (2 * w.ϕ - π))
end

struct MeixnerPollaczek{T,L,P} <: OrthogonalPolynomial{T}
    λ::L
    ϕ::P
end
MeixnerPollaczek{T}(λ::L, ϕ::P) where {T,L,P}=MeixnerPollaczek{T,L,P}(λ,ϕ)
MeixnerPollaczek(λ::L, ϕ::P) where {L,P} = MeixnerPollaczek{promote_type(L, P)}(λ, ϕ)

axes(::MeixnerPollaczek{T}) where {T} = (Inclusion{T}(ℝ), oneto(∞))
orthogonalityweight(S::MeixnerPollaczek{T}) where {T} = MeixnerPollaczekWeight{T}(S.λ, S.ϕ)

jacobimatrix(S::MeixnerPollaczek) = Tridiagonal(RationalVal{1}{2}()*csc(S.ϕ)*oneto(∞),-cot(S.ϕ)*(S.λ:∞), csc(S.ϕ)*(S.λ:∞))

end
