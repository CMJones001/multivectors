module Bivector

import Base: +, -, ==, ^, *, ∘

export V, vec, ∧, e1, e2, e12

struct V{T<:Real} <: Number
    s::T
    e1::T
    e2::T
    e12::T
end

const e1  = V(0, 1, 0, 0)
const e2  = V(0, 0, 1, 0)
const e12 = V(0, 0, 0, 1)

V(s::Real, e1::Real, e2::Real, e12::Real) = V(promote(s,e1,e2,e12)...)
V(a::Real, b::Real) = V(zero(a), a, b, zero(a))
V(a::Real) = V(a, zero(a), zero(a), zero(a))

# Simple operations
+(z::V, w::V) = V(z.s+w.s, z.e1+w.e1, z.e2+w.e2, z.e12+w.e12)
-(z::V, w::V) = V(z.s-w.s, z.e1-w.e1, z.e2-w.e2, z.e12-w.e12)
-(z::V) = V(-z.s, -z.e1, -z.e2, -z.e12)
==(z::V, w::V) = (z.s == w.s) & (z.e1 == w.e1) & (z.e2 == w.e2) & (z.e12 == w.e12) 

# Outer product
∧(z::V, w::V) = V(z.s*w.s,
                  z.s*w.e1 + w.s*z.e1, 
                  z.s*w.e2 + w.s*z.e2, 
                  z.e1*w.e2-z.e2*w.e1 + z.s*w.e12 + w.s*z.e12)

# Geometric Product
# ∘(z::V, w::V) = V(z.e1*w.e1 + z.e2*w.e2 - z.e12*w.e12,
#                   - w.e12*z.e2 + z.e12*w.e2,
#                   + w.e12*z.e1 - z.e12*w.e1,
#                   0)
# *(z::V, w::V) = z∘w + z∧w
*(a::V, b::V) = V(
    a.s*b.s   +  a.e1*b.e1 +  a.e2*b.e2  - a.e12*b.e12,
    a.s*b.e1  +  a.e1*b.s  + a.e12*b.e2  -  a.e2*b.e12,
    a.s*b.e2  +  a.e2*b.s  +  a.e1*b.e12 - a.e12*b.e1,
    a.s*b.e12 + a.e12*b.s  +  a.e1*b.e2  -  a.e2*b.e1
)

Base.:^(z::V, p::Integer) = Base.power_by_squaring(z, p)

# Scalar operations
+(a::Real, z::V) = V(a+z.s, z.e1, z.e2, z.e12)
-(a::Real, z::V) = V(a-z.s, -z.e1, -z.e2, -z.e12)
*(a::Real, z::V) = V(a*z.s, a*z.e1, a*z.e2, a*z.e12)

# We can also give this algebra on 2x2 matricies
E1 = [ 0 +1; +1  0]
E2 = [-1  0;  0 -1]
end
