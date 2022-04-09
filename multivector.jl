module MultivectorComps

import Base: +, -, *, /, ==, show, zero, iszero
export Sca, Vec, Biv, Tri, ∧, ZeroVec

# Classes for the geometric algebra of three dimensions
abstract type Geometric end
struct ZeroVec end
const zv = ZeroVec()

# Vector
struct Vec{T<:Real} <: Geometric
    e1::T
    e2::T
    e3::T
end

# Bivector
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
const v1 = Vec(1, 0, 0)
const v2 = Vec(0, 1, 0)
const v3 = Vec(0, 0, 1)
const v12 = Biv(1, 0, 0)
const v23 = Biv(0, 1, 0)
const v31 = Biv(0, 0, 1)
const I = Tri(1)

const V = Vec
const B = Biv
const T = Tri

V(e1::Real, e2::Real, e3::Real) = V(promote(e1,e2,e3)...)
B(e12::Real, e23::Real, e31::Real) = B(promote(e12,e23,e31)...)

iszero(v::Vec) = iszero(v.e1) & iszero(v.e2) & iszero(v.e3)
iszero(b::Biv) = iszero(b.e12) & iszero(b.e23) & iszero(b.e31)
iszero(t::Tri) = iszero(t.i)
iszero(::ZeroVec) = true

zero(::Vec) = Vec(0, 0, 0)
zero(::Biv) = Biv(0, 0, 0)
zero(::Tri) = Tri(0)

# Generic options (are these slower?)
comp(r::Real) = (r,)
comp(v::Vec) = (v.e1, v.e2, v.e3)
comp(b::Biv) = (b.e12, b.e23, b.e32)
comp(t::Tri) = (t.i,)

add(z::U, w::U) where U<:Geometric = U((comp(z) .+ comp(w))...)
sub(z::U, w::U) where U<:Geometric = U((comp(z) .- comp(w))...)
eq(z::U, w::U) where U<:Geometric = comp(z) == comp(w)

# addition and subtraction
+(z::Vec, w::Vec) = V(z.e1 + w.e1, z.e2 + w.e2, z.e3 + w.e3)
-(z::Vec, w::Vec) = V(z.e1 - w.e1, z.e2 - w.e2, z.e3 - w.e3)
-(z::Vec) = V(-z.e1, -z.e2, -z.e3)
==(z::Vec, w::Vec) = (z.e1 == w.e1) & (z.e2 == w.e2) & (z.e3 == w.e3) 

+(z::Biv, w::Biv) = B(z.e12 + w.e12, z.e23 + w.e23, z.e31 + w.e31)
-(z::Biv, w::Biv) = B(z.e12 - w.e12, z.e23 - w.e23, z.e31 - w.e31)
-(z::Biv) = B(-z.e12, -z.e23, -z.e31)
==(z::Biv, w::Biv) = (z.e12 == w.e12) & (z.e23 == w.e23) & (z.e31 == w.e31) 

+(z::Tri, w::Tri) = Tri(z.i + w.i)
-(z::Tri, w::Tri) = Tri(z.i - w.i)
-(z::Tri) = Tri(-z.i)
==(z::T, w::T) = z.i == w.i

# Zero is the identity under addition
+(z::Geometric, ::ZeroVec) = z
+(::ZeroVec, z::Geometric) = z
-(z::Geometric, ::ZeroVec) = z
-(::ZeroVec, z::Geometric) =-z
+(::ZeroVec, ::ZeroVec) = zv
-(::ZeroVec, ::ZeroVec) = zv

# Zero vector is equal to zero...
==(z::Geometric, ::ZeroVec) = iszero(z)
==(::ZeroVec, z::Geometric) = iszero(z)

# Scalar operations
*(a::Real, z::Vec) = V(z.e1 * a, z.e2 * a, z.e3 * a)
*(z::Vec, a::Real) = V(z.e1 * a, z.e2 * a, z.e3 * a)
/(z::Vec, a::Real) = V(z.e1 / a, z.e2 / a, z.e3 / a)

*(a::Real, z::Biv) = B(z.e12 * a, z.e23 * a, z.e31 * a)
*(z::Biv, a::Real) = B(z.e12 * a, z.e23 * a, z.e31 * a)
/(z::Biv, a::Real) = B(z.e12 / a, z.e23 / a, z.e31 / a)

*(a::Real, z::Tri) = T(z.i * a)
*(z::Tri, a::Real) = T(z.i * a)
/(z::Tri, a::Real) = T(z.i / a)

# Under multiplication, zero vector always returns zero
*(a::Geometric, ::ZeroVec) = zv
*(z::ZeroVec, a::Geometric) = zv
*(::ZeroVec, ::ZeroVec) = zv
*(::ZeroVec, r::Real) = zv
*(r::Real, z::ZeroVec) = zv

# Misc options
show(io::IO, v::Vec) = print(io, "$(v.e1)×e1 + $(v.e2)×e2 + $(v.e3)×e3")
show(io::IO, v::Biv) = print(io, "$(v.e12)×e1∧e2 + $(v.e23)×e2∧e3 + $(v.e31)×e3∧e1")
show(io::IO, v::Tri) = print(io, "$(v.i)×e1∧e2∧e3")
show(io::IO, ::ZeroVec) = print(io, "∅")

# Outer product with scalars is just multiplication
∧(z::Geometric, a::Real) = z*a
∧(a::Real, z::Geometric) = z*a

# Outer product with Nothing gives nothing
∧(z::Geometric, ::ZeroVec) = zv
∧(::ZeroVec, ::Geometric) = zv

# Outer product of vectors gives a bivector
∧(z::Vec, w::Vec) = Biv(
    z.e1 * w.e2 - z.e2 * w.e1,
    z.e2 * w.e3 - z.e3 * w.e2,
    z.e3 * w.e1 - z.e1 * w.e3,
)

∧(b::Biv, v::Vec) = Tri(b.e12 * v.e3 + b.e23 * v.e1 + b.e31 * v.e2)
∧(v::Vec, b::Biv) = Tri(b.e12 * v.e3 + b.e23 * v.e1 + b.e31 * v.e2)
∧(::ZeroVec, ::ZeroVec) = zv

end

module Multivector
using ..MultivectorComps

import Base: +, -, *, /, ==, show, zero, iszero
import ..MultivectorComps: V, B, T, v1, v2, v3, v12, v23, v31, I, ∧, zv, Geometric
export V, B, T, v1, v2, v3, v12, v23, v31, I
export Vec, Biv, Tri, ∧, Mul, zv

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
Mul() = Mul{Real}(zv, zv, zv, zv)
# We need to define a generic type when all values are zero still
Mul(::ZeroVec, ::ZeroVec, ::ZeroVec, ::ZeroVec) = Mul{UInt8}(zv, zv, zv, zv)

==(m::Mul, s::Real) = (m.s == s) & iszero(m.b) & iszero(m.v) & iszero(m.p)
==(b::Real, m::Mul) = (m.s == s) & iszero(m.b) & iszero(m.v) & iszero(m.p)
==(m::Mul, v::Vec) = (m.v == v) & iszero(m.s) & iszero(m.b) & iszero(m.p)
==(v::Vec, m::Mul) = (m.v == v) & iszero(m.s) & iszero(m.b) & iszero(m.p)
==(m::Mul, b::Biv) = (m.b == b) & iszero(m.s) & iszero(m.v) & iszero(m.p)
==(b::Biv, m::Mul) = (m.b == b) & iszero(m.s) & iszero(m.v) & iszero(m.p)
==(m::Mul, p::Tri) = (m.p == p) & iszero(m.b) & iszero(m.v) & iszero(m.s)
==(p::Tri, m::Mul) = (m.p == p) & iszero(m.b) & iszero(m.v) & iszero(m.s)

# Convert to multivectors if no other addition rules are available
+(z::Geometric, w::Geometric) = Mul(z) + Mul(w)
-(z::Geometric, w::Geometric) = Mul(z) - Mul(w)

# This is likely to cause problems elsewhere?!
# Maybe replace with a NoVec type
+(z::Mul, w::Mul) = Mul(z.s + w.s, z.v + w.v, z.b + w.b, z.p + w.p)
-(z::Mul, w::Mul) = Mul(z.s - w.s, z.v - w.v, z.b - w.b, z.p - w.p)
==(z::Mul, w::Mul) = (z.s == w.s) & (z.v == w.v) & (z.b == w.b) & (z.p == w.p)
iszero(m::Mul) = iszero(m.s) & iszero(m.v) & iszero(m.b) & iszero(m.p)

∧(z::Mul, w::Mul) = Mul(
    z.s * w.s,
    z.s * w.v + z.v * w.s,
    z.s * w.b + z.b * z.s + z.v ∧ w.v,
    z.s * w.p + z.p * w.s + z.v ∧ w.b + z.b ∧ w.v
)

end
