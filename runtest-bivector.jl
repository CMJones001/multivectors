module test

include("bivectors.jl")
using .Bivector
using Test

@testset "Outer product" begin
    @testset "Orthogonal Vectors" begin
        a = V(1, 0)
        b = V(0, 1)
        @test a∧b == V(0, 0, 0, 1)

        a = e1
        b = e2
        @test a∧b == e12

        # Test with mixed values
        a = V(1, 0.0)
        b = V(0, 1)
        @test a∧b == V(0, 0, 0, 1)

        a = V(2, 0)
        b = V(0, 3)
        @test a∧b == V(0, 0, 0, 6)
    end

    @testset "Parallel vectors" begin
        test_zero_outer(a, b) = @test (a ∧ b) == V(0)

        test_zero_outer(V(1, 0), V(1, 0))
        test_zero_outer(V(2, 0), V(1, 0))
        test_zero_outer(V(1, 1), V(1, 1))
        test_zero_outer(V(1, 1), V(3, 3))
    end

    @testset "Arbitary vectors" begin
        @test V(2, 3) ∧ V(1,-1) == -5*e12
        @test V(4,-2) ∧ V(3, 4) == (16+6)*e12
    end

    @testset "Scalar product" begin
        @test V(1, 2) ∧ V(3) == V(3, 6)
        @test V(4) ∧ V(2, 3) == V(8, 12)
    end

    @testset "With bivectors" begin
        # Any vector product with a bivector should be zero
        @test e12 ∧ e1 == V(0)
        @test e12 ∧ (e1 + e2) == V(0)
        @test 2*e12 ∧ V(2, 4) == V(0)

        # However, I think scalar product should still work?
        @test e12 ∧ V(2) == 2*e12
    end
end

@testset "Geometric product" begin
    @testset "Orthogonal => outer product" begin
        a = e1
        b = e2
        @test a * b == a ∧ b
        @test a * b == e12

        a = e1
        b = 2*e2
        @test a * b == a ∧ b
        @test a * b == 2*e12
    end

    @testset "Parallel vectors => dot product" begin
        @test e1 * e1 == V(1)
        @test (4*e1) * (2*e1) == V(8)
    end

    @testset "Arbitary vectors" begin
        a = V(1, 2)
        b = V(3, 4)
        @test a * b == V(11, 0, 0, -2)
        @test a * b == 11 - 2*e12
    end
end

@testset "Rotations" begin
    @test e12*e1 == -e2
    @test e12*e2 == +e1
    @test e1*e12 == +e2
    @test e2*e12 == -e1

    @test e12*e12 == V(-1)
    @test e12^2 == V(-1)
end

end
