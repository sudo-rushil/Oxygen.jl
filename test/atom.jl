using Oxygen

println("Testing Atom")

@testset "OxygenAtom" begin
    carbon = OxygenAtom(6, "C", 12.011, Oxygen.Electrons(6), 2.55)
    @test typeof(carbon.number) == Int
    @test typeof(carbon.symbol) == String
    @test typeof(carbon.mass) == Float64
    @test typeof(carbon.orbitals) == Oxygen.Electrons
    @test typeof(carbon.electronegativity) == Float64
    @test carbon == OxygenAtom(6)
    @test carbon == OxygenAtom("C")
end

@testset "Printing" begin
    carbon = OxygenAtom(6, "C", 12.011, Oxygen.Electrons(6), 2.55)
    @test repr(carbon) == "\"C\" atom (number 6)"
end
