module MultivectorComps

import Base: +, -, *, /, ==, show, zero, iszero, |
export Sca, Vec, Biv, Tri, ∧, ZeroVec, norm, contract

"""
Core components required for geometric algebra in three dimensions, ``G_3``.
"""

abstract type Geometric end
struct ZeroVec end
const zv = ZeroVec()

"""
Vectors, or 1 dimensional elements.

Given in terms of three orthogonal basis vectors ``{e1, e2, e3}``.
"""
struct Vec{T<:Real} <: Geometric
    e1::T
    e2::T
    e3::T
end

"""
Bivectors, or 2 dimensional elements

Given in terms of the basis 2-blades ``{e1 ∧ e2, e2 ∧ e3, e3 ∧ e1}``.
"""
struct Biv{T<:Real} <: Geometric
    e12::T
    e23::T
    e31::T
end

# Pseudoscalar (trivector) element
struct Tri{T<:Real} <: Geometric
    i::T
end

# Use v as the notation for the basis vectors, as e is used for exponentiation
const v1  = Vec(1, 0, 0)
const v2  = Vec(0, 1, 0)
const v3  = Vec(0, 0, 1)
const v12 = Biv(1, 0, 0)
const v23 = Biv(0, 1, 0)
const v31 = Biv(0, 0, 1)
const I   = Tri(1)

# Subscripts are allowed, but they are awkward to type, and don't work for bivectors
const e₁ = Vec(1, 0, 0)
const e₂ = Vec(0, 1, 0)
const e₃ = Vec(0, 0, 1)

const V = Vec
const B = Biv
const T = Tri

V(e1::Real, e2::Real, e3::Real) = V(promote(e1,e2,e3)...)
B(e12::Real, e23::Real, e31::Real) = B(promote(e12,e23,e31)...)

# Generic options, compiler makes this as fast as handwritten version 
comp(r::Real) = (r,)
comp(v::Vec)  = (v.e1, v.e2, v.e3)
comp(b::Biv)  = (b.e12, b.e23, b.e31)
comp(t::Tri)  = (t.i,)


# Simple scalar operations ####################################################

# Addition, subtraction, and negation between like types
+(z::U, w::U) where U<:Geometric = U((comp(z) .+ comp(w))...)
-(z::U, w::U) where U<:Geometric = U((comp(z) .- comp(w))...)
-(z::U) where U<:Geometric = U(map(-, comp(z))...)
==(z::U, w::U) where U<:Geometric = comp(z) == comp(w)

# Scalar multiplication
*(a::Real, z::U) where U<:Geometric = U((a .* comp(z))...)
*(z::U, a::Real) where U<:Geometric = U((a .* comp(z))...)
/(z::U, a::Real) where U<:Geometric = U((a ./ comp(z))...)

iszero(v::Geometric) = all(iszero.(comp(v)))
iszero(::ZeroVec) = true

zero(::Vec) = Vec(0, 0, 0)
zero(::Biv) = Biv(0, 0, 0)
zero(::Tri) = Tri(0)

# Zero is the identity under addition
+(z::Geometric, ::ZeroVec) = z
+(::ZeroVec, z::Geometric) = z
-(z::Geometric, ::ZeroVec) = z
-(::ZeroVec, z::Geometric) = -z
+(::ZeroVec, ::ZeroVec)    = zv
-(::ZeroVec, ::ZeroVec)    = zv

# Zero vector is equal to zero...
==(z::Geometric, ::ZeroVec) = iszero(z)
==(::ZeroVec, z::Geometric) = iszero(z)

# Under multiplication, zero vector always returns zero
*(a::Geometric, ::ZeroVec)  = zv
*(z::ZeroVec, a::Geometric) = zv
*(::ZeroVec, ::ZeroVec)     = zv
*(::ZeroVec, r::Real)       = zv
*(r::Real, z::ZeroVec)      = zv

# Misc options
show(io::IO, v::Vec)    = print(io, "$(v.e1)×e₁ + $(v.e2)×e₂ + $(v.e3)×e₃")
show(io::IO, v::Biv)    = print(io, "$(v.e12)×e₁∧e₂ + $(v.e23)×e₂∧e₃ + $(v.e31)×e₃∧e₁")
show(io::IO, v::Tri)    = print(io, "$(v.i)×e₁∧e₂∧e₃")
show(io::IO, ::ZeroVec) = print(io, "∅")

# Outer products ##############################################################

# Outer product with scalars is just multiplication
∧(z::Geometric, a::Real) = z*a
∧(a::Real, z::Geometric) = z*a

# Outer product with Nothing gives nothing
∧(z::Geometric, ::ZeroVec) = zv
∧(::ZeroVec, ::Geometric)  = zv

# Outer product of vectors gives a bivector
∧(z::Vec, w::Vec) = Biv(
    z.e1 * w.e2 - z.e2 * w.e1,
    z.e2 * w.e3 - z.e3 * w.e2,
    z.e3 * w.e1 - z.e1 * w.e3,
)

∧(b::Biv, v::Vec) = Tri(b.e12 * v.e3 + b.e23 * v.e1 + b.e31 * v.e2)
∧(v::Vec, b::Biv) = Tri(b.e12 * v.e3 + b.e23 * v.e1 + b.e31 * v.e2)
∧(::ZeroVec, ::ZeroVec) = zv


# Contractions ################################################################

# Dot producs
|(z::Vec, w::Vec) = sum(comp(z) .* comp(w))
|(z::Real, w::Real) = z*w

norm(z::Geometric) = √(z|z)

# Trivial cases
contract(a::Real, b::Real) = a * b
contract(a::Biv, b::Biv)   = -(a.e12 * b.e12 + a.e23 * b.e23 + a.e31 * b.e31)
contract(a::Tri, b::Tri)   = -a.i * b.i

contract(a::Real, b::Geometric) = a * b
contract(a::Vec, b::Vec) = a | b

# All cases of grade(a) > grade(b) = 0
contract(a::Vec, b::Real) = zv
contract(a::Biv, b::Real) = zv
contract(a::Biv, b::Vec)  = zv
contract(a::Tri, b::Real) = zv
contract(a::Tri, b::Vec)  = zv
contract(a::Tri, b::Biv)  = zv

function contract(a::Vec, b::Biv)
    Vec(
        b.e31 * a.e3 - b.e12 * a.e2,
        b.e12 * a.e1 - b.e23 * a.e3,
        b.e23 * a.e2 - b.e31 * a.e1
    )
end

contract(a::Vec, b::Tri) = -b.i * Biv(a.e3, a.e1, a.e2)
contract(a::Biv, b::Tri) = -b.i * Vec(a.e23, a.e31, a.e12)

const << = contract

end

module Multivectors
using ..MultivectorComps

import Base: +, -, *, /, ==, show, zero, iszero
import ..MultivectorComps: V, B, T, v1, v2, v3, v12, v23, v31, I, ∧, zv, Geometric
export V, B, T, v1, v2, v3, v12, v23, v31, I, norm
export Vec, Biv, Tri, ∧, Mul, zv, contract

struct Mul{T<:Real}
    s::Union{T, ZeroVec}
    v::Union{Vec{T}, ZeroVec}
    b::Union{Biv{T}, ZeroVec}
    p::Union{Tri{T}, ZeroVec}
end

∧(z::Mul, w::Geometric) = z ∧ Mul(w)
∧(z::Geometric, w::Mul) = Mul(z) ∧ w

Mul(s::Real) = Mul(s, zv, zv, zv)
Mul(v::Vec)  = Mul(zv, v, zv, zv)
Mul(b::Biv)  = Mul(zv, zv, b, zv)
Mul(t::Tri)  = Mul(zv, zv, zv, t)
# We need to define a generic type when all values are zero still
Mul() = Mul{UInt8}(zv, zv, zv, zv)
Mul(::ZeroVec, ::ZeroVec, ::ZeroVec, ::ZeroVec) = Mul{UInt8}(zv, zv, zv, zv)

==(m::Mul, s::Real) = (m.s == s) & iszero(m.b) & iszero(m.v) & iszero(m.p)
==(b::Real, m::Mul) = (m.s == s) & iszero(m.b) & iszero(m.v) & iszero(m.p)
==(m::Mul, v::Vec)  = (m.v == v) & iszero(m.s) & iszero(m.b) & iszero(m.p)
==(v::Vec, m::Mul)  = (m.v == v) & iszero(m.s) & iszero(m.b) & iszero(m.p)
==(m::Mul, b::Biv)  = (m.b == b) & iszero(m.s) & iszero(m.v) & iszero(m.p)
==(b::Biv, m::Mul)  = (m.b == b) & iszero(m.s) & iszero(m.v) & iszero(m.p)
==(m::Mul, p::Tri)  = (m.p == p) & iszero(m.b) & iszero(m.v) & iszero(m.s)
==(p::Tri, m::Mul)  = (m.p == p) & iszero(m.b) & iszero(m.v) & iszero(m.s)

# Convert to multivectors if no other addition rules are available
+(z::Geometric, w::Geometric) = Mul(z) + Mul(w)
-(z::Geometric, w::Geometric) = Mul(z) - Mul(w)

# This is likely to cause problems elsewhere?!
# Maybe replace with a NoVec type
comp(m::Mul) = (m.s, m.v, m.b, m.p)
+(z::Mul, w::Mul)  = Mul((comp(z) .+ comp(w))...)
-(z::Mul, w::Mul)  = Mul((comp(z) .- comp(w))...)
==(z::Mul, w::Mul) = comp(z) == comp(w)
iszero(m::Mul)     = all(iszero.(comp(m)))

∧(z::Mul, w::Mul) = Mul(
    z.s * w.s,
    z.s * w.v + z.v * w.s,
    z.s * w.b + z.b * z.s + z.v ∧ w.v,
    z.s * w.p + z.p * w.s + z.v ∧ w.b + z.b ∧ w.v
)

end
