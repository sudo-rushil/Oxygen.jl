using Oxygen

println("Testing Molecules")

@testset "OxygenMol" begin
    water = OxygenMol([OxygenAtom("H"), OxygenAtom("O"), OxygenAtom("H")], [[(2, 1)], [(1, 1), (3, 1)], [(3, 1)]])
    @test typeof(water.atoms) == Array{OxygenAtom,1}
    @test typeof(water.adj) <: Array
    morewater = OxygenMol([1, 8, 1])
    morewater.adj = [[(2, 1)], [(1, 1), (3, 1)], [(3, 1)]]
    @test water == morewater
    alsowater = OxygenMol(["H", "O", "H"])
    alsowater.adj = [[(2, 1)], [(1, 1), (3, 1)], [(3, 1)]]
    @test water == alsowater
end
