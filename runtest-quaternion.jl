module test
include("quaternion.jl")

using .Quat
using Test

@testset "Test fundamental operations" begin
   @test ii * ii ≈ -rr
   @test jj * jj ≈ -rr
   @test kk * kk ≈ -rr
   @test ii * jj * kk ≈ -rr
end

@testset "Test identities" begin
    @test jj * kk ≈ +ii
    @test kk * jj ≈ -ii

    @test kk * ii ≈ +jj
    @test ii * kk ≈ -jj

    @test ii * jj ≈ +kk
    @test jj * ii ≈ -kk

    @test ii * rr ≈ ii
    @test rr * ii ≈ ii

    @test jj * rr ≈ jj
    @test rr * jj ≈ jj

    @test kk * rr ≈ kk
    @test rr * kk ≈ kk
end

@testset "Scalar Operations" begin
    @test 3*Q(1, 2, 3, 4) ≈ Q(3, 6, 9, 12)
    @test Q(1, 2, 3, 4)*3.0 ≈ Q(3, 6, 9, 12)
    @test Q(1, 2, 3, 4) ≈ Q(10, 20, 30, 40)/10.0

    @testset "Equality" begin
        @test Q(1, 2, 3, 4) ≈ 1 + 2ii + 3jj + 4kk
        @test Q(1, 0, 0, 0) ≈ 1
        @test !(Q(1, 0, 2, 0) ≈ 1)
    end

    @testset "Alternate notations" begin
        @test 3.0*Q(1, 2, 3, 4) ≈ Q(3, 6, 9, 12)
        @test 3*Q(1, 2, 3, 4) ≈ Q(3, 6, 9, 12)
        @test 3Q(1, 2, 3, 4) ≈ Q(3, 6, 9, 12)
    end
end


@testset "Test non-trivial multiplication" begin
    @testset "interseting results" begin
        @test (ii + jj)*(ii + jj) ≈ -2.0rr
    end

    @testset "Arbitary values" begin
        z = Q(1, 2, 3, 4)
        w = Q(5, 6, 7, 8)
        @test z * w ≈ Q(-60, 12, 30, 24)
        @test w * z ≈ Q(-60, 20, 14, 32)

        z = Q(+8, +4, -2, +6)
        w = Q(-1, -4, -4, +2)
        @test z * w ≈ Q(-12, -16, -62, -14)
        @test w * z ≈ Q(-12, -56,  +2, +34)
    end
end

@testset "Powers" begin
    @testset "Identities" begin
        @test ii^2 ≈ -rr
        @test jj^2 ≈ -rr
        @test kk^2 ≈ -rr
    end

    @testset "Complex squaring" begin
        @test (ii+jj)^2 ≈ -2rr
        @test (ii+jj+kk)^2 ≈ -3rr
        @test Q(1, 0, 3, -1)^2 ≈ Q(-9, 0, 6, -2)
        @test Q(1, 0, 3, -1)^2 ≈ -9 + 6jj - 2kk
    end
    
    @testset "Higher powers" begin
        @test Q(1, 0, 3, -1)^3 ≈ Q(-29, 0, -21, +7)
        @test Q(1, 0, 3, -1)^3 ≈ -29rr - 21jj + 7kk
    end
end

end
