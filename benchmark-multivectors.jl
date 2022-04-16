using BenchmarkTools
include("multivector.jl")
using .Multivectors
using Printf

"""
Benchmarking comparison of the specific vector operations and their
multivector representations.

We see that the specific operations are far faster, and should be
used as much as possible.
"""

N = 2^20

# Generate random Vectors and their multivector representation
new_vectors(n::Integer) = [V(rand(-10:10, 3)...) for i in 1:n]
new_multivectors(n) = map(Mul, new_vectors(n))

# Print the average time taken for each outer product 
function print_time(results, label, n::Integer)
    mean_time = median(results.times) / n
    @printf("Mean time for %s is %.2f ns (N=%d).\n", label, mean_time, n)
end

P = new_vectors(N)
Q = new_vectors(N)

res = @benchmark $P .∧ $Q
print_time(res, "simple", N)

Pm = new_multivectors(N)
Qm = new_multivectors(N)

res = @benchmark $Pm .∧ $Qm
print_time(res, "multivector", N)
