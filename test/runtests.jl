using Ellasy
using Test
using Aqua

@testset "Ellasy.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(Ellasy)
    end
    # Write your tests here.
end
