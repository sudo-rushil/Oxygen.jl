using Oxygen, Test

@testset "All Tests" begin
    include("periodic.jl")
    include("electrons.jl")
    include("atom.jl")
end
