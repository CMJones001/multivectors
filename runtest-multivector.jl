module test
include("multivector.jl")

using .Multivector
using Test
import .Multivector: V, B, T, v1, v2, v3, v12, v23, v31, I

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
end
