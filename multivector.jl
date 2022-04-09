module Multivector

import Base: +, -, *, /, ==, show, zero, iszero
export Vec, Biv, Tri, ∧

# Classes for the geometric algebra of three dimensions

# Vector
struct Vec{T<:Real} <: Number
    e1::T
    e2::T
    e3::T
end

# Bivector
struct Biv{T<:Real} <: Number
    e12::T
    e23::T
    e31::T
end

# Pseudoscalar (trivector) element
struct Tri{T<:Real}
    i::T
end

# Combination of all graded elements
struct Mul{T<:Real} 
    s::T
    v::Vec{T}
    b::Biv{T}
    p::T
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
# T(i::Real) = T(promote(i)...)
iszero(v::Vec) = iszero(v.e1) & iszero(v.e2) & iszero(v.e3)
iszero(b::Biv) = iszero(b.e12) & iszero(b.e23) & iszero(b.e31)
iszero(t::Tri) = iszero(t.i)


# Simple operations
+(z::Vec, w::Vec) = V(z.e1 + w.e1, z.e2 + w.e2, z.e3 + w.e3)
-(z::Vec, w::Vec) = V(z.e1 - w.e1, z.e2 - w.e2, z.e3 - w.e3)
-(z::Vec) = V(-z.e1, -z.e2, -z.e3)
==(z::Vec, w::Vec) = (z.e1 == w.e1) & (z.e2 == w.e2) & (z.e3 == w.e3) 

+(z::Biv, w::Biv) = B(z.e12 + w.e12, z.e23 + w.e23, z.e31 + w.e31)
-(z::Biv, w::Biv) = B(z.e12 - w.e12, z.e23 - w.e23, z.e31 - w.e31)
-(z::Biv) = B(-z.e12, -z.e23, -z.e31)
==(z::Biv, w::Biv) = (z.e12 == w.e12) & (z.e23 == w.e23) & (z.e31 == w.e31) 

+(z::Tri, w::Tri) = T(z.i + w.i)
-(z::Tri, w::Tri) = T(z.i - w.i)
-(z::Tri) = T(-z.i)

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

# Misc options
show(io::IO, v::Vec) = print(io, "$(v.e1)×e1 + $(v.e2)×e2 + $(v.e3)×e3")
show(io::IO, v::Biv) = print(io, "$(v.e12)×e1∧e2 + $(v.e23)×e2∧e3 + $(v.e31)×e3∧e1")
show(io::IO, v::Tri) = print(io, "$(v.i)×e1∧e2∧e3")

# Outer product with scalars is just multiplication
∧(z::Vec, a::Real) = z*a
∧(a::Real, z::Vec) = z*a

∧(z::Biv, a::Real) = z*a
∧(a::Real, z::Biv) = z*a

∧(z::Tri, a::Real) = z*a
∧(a::Real, z::Tri) = z*a

# Outer product of vectors gives a bivector
∧(z::Vec, w::Vec) = Biv(
    z.e1 * w.e2 - z.e2 * w.e1,
    z.e2 * w.e3 - z.e3 * w.e2,
    z.e3 * w.e1 - z.e1 * w.e3,
)

∧(b::Biv, v::Vec) = Tri(b.e12 * v.e3 + b.e23 * v.e1 + b.e31 * v.e2)
∧(v::Vec, b::Biv) = Tri(b.e12 * v.e3 + b.e23 * v.e1 + b.e31 * v.e2)

end
