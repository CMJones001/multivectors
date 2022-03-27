module Quat

export Quaternion, rr, ii, jj, kk, Q

struct Quaternion<: Number
    r::Float64
    i::Float64
    j::Float64
    k::Float64
end

Q(r, i, j, k) = Quaternion(r, i, j, k)
const rr = Quaternion(1.0, 0.0, 0.0, 0.0)
const ii = Quaternion(0.0, 1.0, 0.0, 0.0)
const jj = Quaternion(0.0, 0.0, 1.0, 0.0)
const kk = Quaternion(0.0, 0.0, 0.0, 1.0)

# Equality
Base.:≈(z::Quaternion, w::Quaternion) = z.r ≈ w.r && z.i ≈ w.i && z.j ≈ w.j && z.k ≈ w.k
Base.:≈(z::Quaternion, a::Real) = z.r ≈ a && z.i ≈ 0.0 && z.j ≈ 0.0 && z.k ≈ 0.0
Base.:≈(a::Real, z::Quaternion) = z.r ≈ a && z.i ≈ 0.0 && z.j ≈ 0.0 && z.k ≈ 0.0
Base.:isequal(z::Quaternion, w::Quaternion) = isequal(z.r, w.r) && isequal(z.i, w.i) && isequal(z.j, w.j) && isequal(z.k, w.k)

# Quaternion operations
Base.:+(z::Quaternion, w::Quaternion) = Quaternion(z.r + w.r, z.i + w.i, z.j + w.j, z.k + w.k) 
Base.:-(z::Quaternion, w::Quaternion) = Quaternion(z.r - w.r, z.i - w.i, z.j - w.j, z.k - w.k) 
Base.:*(z::Quaternion, w::Quaternion) = Quaternion(
    z.r*w.r - z.i*w.i - z.j*w.j - z.k*w.k,
    z.r*w.i + z.i*w.r + z.j*w.k - z.k*w.j,
    z.r*w.j - z.i*w.k + z.j*w.r + z.k*w.i,
    z.r*w.k + z.i*w.j - z.j*w.i + z.k*w.r,
)
Base.:^(z::Quaternion, p::Integer) = Base.power_by_squaring(z, p)

# Negation
Base.:-(z::Quaternion) = Quaternion(-z.r, -z.i, -z.j, -z.k) 

# Scalar multiplication
Base.:*(a::Real, z::Quaternion) = Quaternion(a*z.r, a*z.i, a*z.j, a*z.k)
Base.:*(z::Quaternion, a::Real) = Quaternion(a*z.r, a*z.i, a*z.j, a*z.k)
Base.:/(z::Quaternion, a::Float64) = Quaternion(z.r/a, z.i/a, z.j/a, z.k/a)

Base.:+(a::Real, z::Quaternion) = Quaternion(a+z.r, z.i, z.j, z.k)
Base.:-(a::Real, z::Quaternion) = Quaternion(a-z.r, z.i, z.j, z.k)
Base.:+(z::Quaternion, a::Real) = Quaternion(a+z.r, z.i, z.j, z.k)
Base.:-(z::Quaternion, a::Real) = Quaternion(z.r-a, z.i, z.j, z.k)

z = Quaternion(1, 2, 3, 4)
w = Quaternion(5, 6, 7, 8)

end
