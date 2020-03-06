using Oxygen

println("Testing Atom")

@testset "StaticAtom" begin
    carbon = StaticAtom(6, "C", 12.011, 2.55)
    @test typeof(carbon.number) == Int
    @test typeof(carbon.symbol) == String
    @test typeof(carbon.mass) == Float64
    @test typeof(carbon.electronegativity) == Float64
    @test carbon == StaticAtom(6)
    @test carbon == StaticAtom("C")
end

@testset "Atom" begin
    carbon = OxygenAtom(6, "C", 12.011, 2.55, 2, 0, 4)
    @test typeof(carbon.number) == Int
    @test typeof(carbon.symbol) == String
    @test typeof(carbon.mass) == Float64
    @test typeof(carbon.electronegativity) == Float64
    @test typeof(carbon.formalcharge) == Int
    @test typeof(carbon.radicals) == Int
    @test typeof(carbon.degree) == Int
    @test Oxygen.atom_equal(carbon, OxygenAtom(6))
    @test Oxygen.atom_equal(carbon, OxygenAtom("C"))
end
