module test
include("multivector.jl")

using .Multivectors
using Test
import .Multivectors: V, B, T, v1, v2, v3, v12, v23, v31, I

@testset "Basic scalar opts" begin
    @testset "Vectors" begin
        z = Vec(2, 3, 4)
        α = 2
        @test Vec(4, 6, 8) == z * α

        # Shorter notation
        z = V(1, 0, 1)
        α = 3
        @test V(3, 0, 3) == z * α
        @test V(3, 0, 3) == z ∧ α

        # Basis vectors
        @test V(4, 2, 3) == 4v1 + V(0, 2, 3)
        @test V(4, 2, 3) == 4v1 + 2v2 + 3v3
    end

    @testset "Bivectors" begin
        z = Biv(2, 3, 4)
        α = 2
        @test Biv(4, 6, 8) == z * α

        # Shorter notation
        z = B(1, 0, 1)
        α = 3
        @test B(3, 0, 3) == z * α
        @test B(3, 0, 3) == z ∧ α

        # Basis vectors
        @test B(4, 2, 3) == 4v12 + B(0, 2, 3)
        @test B(4, 2, 3) == 4v12 + 2v23 + 3v31

        @test 3v1 ∧ 3v2 - 3v12 isa Biv
        @test 3v1 ∧ 3v2 - 3v12 - 3v23 == 3(2v12 - v23)
    end
end

@testset "Outer product of vectors" begin
    @testset "basis vectors" begin
        @test v1 ∧ v2 == v12
        @test iszero(v1 ∧ v1) 
        @test v2 ∧ v1 == -v12
        @test v2 ∧ v3 == v23
    end

    @testset "arbitary vectors" begin
        @test V(1, 0, 2) ∧ V(0, 3, 1) == B(3, -6, -1)
        @test (v1 + 2v3) ∧ (3v2 + v3) == 3v12 - 6v23 - v31
        @test V(1, -1, 2) ∧ V(2, 4, 0) == B(6, -8, 4)
    end

    @testset "basis vectors + bivector" begin
        @test iszero(v1 ∧ v12) 
        @test iszero(v1 ∧ v31) 
        @test v2 ∧ v31 == I
        @test v31 ∧ v2 == I

        @test B(1, 3, 4) ∧ V(1, 0, -1) == T(2)
        @test (1v12 + 3v23 + 4v31) ∧ (v1 - v3) == 2I
    end
end

@testset "Dot product" begin
    @testset "vectors" begin
        @test v1 | v2 == 0
        @test v1 | v1 == 1

        @test V(1, 2, 3) | v2 == 2
        @test V(3, 4, 2) | V(1, 0, 2) == 7
    end
end

@testset "Misc functions" begin
    @testset "vectors" begin
        @test norm(v1) == 1
        @test norm(v2) == 1
        @test norm(V(3, 4, 0)) == 5.0
    end
end

@testset "Contractions" begin
    @testset "Trivial cases" begin
        # Scalars
        @test contract(3, 2) == 6
        @test contract(0, 2) == 0

        # Vectors
        @test contract(v1, v2) == 0
        @test contract(v1, v1) == 1
        @test contract(v1, 5v1) == 5
        @test contract(V(2, 3, 1), V(3, 0, 1)) == 7

        # Mismatched ranks
        @test contract(v1, 3) == zv

        # Scalar and vector comps
        @test contract(3, v1) == 3v1
        @test contract(3, 2v2) == 6v2
        @test contract(4, v31) == 4v31
    end

    @testset "Vector to Bivector" begin
        @test contract(2v1, v12) == +2v2
        @test contract(2v2, v12) == -2v1
        @test contract(3v3+2v2, 5v23) == Vec(0, -15, 10)
    end

    @testset "Vector to Trivector" begin
        @test contract(v23, I) == -1v1
        @test contract(-2v12, 3I) == 6v3
        @test contract(4v31, I) == -4v2
    end

    @testset "Bivector to bivector" begin
        @test contract(v12, v12) == -1
        @test contract(v12, v23) == 0
        @test contract(v23, v31 + 2v23) == -2
    end

    @testset "Bivector to Tri" begin
        @test contract(v12, I) == -v3
        @test contract(3v31, 3I) == -9v2
        @test contract(v31 + 2v12, 6I) == -6v2 - 12v3
    end

    @testset "Trivectors to trivectors" begin
        @test contract(I, 3I) == -3
    end
end

@testset "Multivector General" begin
    @test Mul(2) ∧ Mul(3) == Mul(6)

    @testset "Vector combinations" begin
        @test Mul(2) ∧ Mul(v1) == Mul(2v1)
        @test Mul(v1) ∧ Mul(v2) == Mul(v12)
        @test Mul(4v2) ∧ Mul(v1) == Mul(-4v12)
        @test Mul(v1) + Mul(2v2)  == Mul(V(1, 2, 0))
    end

    @testset "Parallel vectors" begin
        @test iszero(Mul(v2) ∧ Mul(v2)) 
        @test Mul(v2) ∧ Mul(v2) == Mul()
        @test Mul(v1) ∧ Mul(4v1) == Mul()
    end

    @testset "Multivectors" begin
        @test Mul(2) ∧ Mul(4v31) == Mul(8v31)
        @test Mul(v31) ∧ Mul(v2) == Mul(Tri(1))
        @test Mul(v1) ∧ Mul(v23) == Mul(I)
        @test Mul(4v1) ∧ Mul(2v23) == Mul(8I)

        @test Mul(v31) ∧ Mul(v1) == Mul()
        @test Mul(v31) ∧ Mul(32) == Mul()
    end

    @testset "Cross components" begin
        @test Mul(2) ∧ Mul(v1) == 2v1
        @test Mul(4v1) ∧ Mul(2v23) == 8I
        @test Mul(4v1) ∧ 2v23 == 8I
        @test Mul(4v1) ∧ 2v23 == 8I
    end

    @testset "Cross addition" begin
        @test v1 + I == Mul(zv, v1, zv, I)
        @test 3v1 + 4v2 - 2I == Mul(zv, V(3, 4, 0), zv, -2I)
    end
end
end
