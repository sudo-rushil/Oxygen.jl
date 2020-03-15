using Oxygen

println("Testing Electrons")

@testset "Electrons" begin
    carbon = Oxygen.parse_config("[He]2s2 2p2")
    @test typeof(carbon.root) == String
    @test typeof(carbon.orbitals) <: Tuple
    @test carbon == Oxygen.Electrons(6)
    @test carbon == Oxygen.Electrons("C")
end

@testset "Printing" begin
    carbon = Oxygen.Electrons(6)
    @test repr(carbon) == "[He] ((\"2s\", 2), (\"2p\", 2))"
end
