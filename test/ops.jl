using Oxygen

println("Testing Operations")

@testset "MolOps" begin
   ethyl_acetate = smilestomol("CC(=O)OCC")
   @test Oxygen.getbonds(ethyl_acetate) == [
      (1, 2, 1.0)
      (2, 3, 2.0)
      (2, 4, 1.0)
      (4, 5, 1.0)
      (5, 6, 1.0)
   ]
   @test Oxygen.getadjmatrix(ethyl_acetate) == [
      0.0 1.0 0.0 0.0 0.0 0.0
      1.0 0.0 2.0 1.0 0.0 0.0
      0.0 2.0 0.0 0.0 0.0 0.0
      0.0 1.0 0.0 0.0 1.0 0.0
      0.0 0.0 0.0 1.0 0.0 1.0
      0.0 0.0 0.0 0.0 1.0 0.0
   ]
   @test Oxygen.valencies(ethyl_acetate) == [1, 3, 1, 2, 2, 1]
   @test Oxygen.atomtypes(ethyl_acetate) == [
      1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
      1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
      0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
      0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
      1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
      1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
   ]

   benzene = smilestomol("c1ccccc1")
   @test Oxygen.getbonds(benzene) == [
      (1, 2, 1.5)
      (1, 6, 1.5)
      (2, 3, 1.5)
      (3, 4, 1.5)
      (4, 5, 1.5)
      (5, 6, 1.5)
   ]
   @test Oxygen.getadjmatrix(benzene) == [
      0.0  1.5  0.0  0.0  0.0  1.5
      1.5  0.0  1.5  0.0  0.0  0.0
      0.0  1.5  0.0  1.5  0.0  0.0
      0.0  0.0  1.5  0.0  1.5  0.0
      0.0  0.0  0.0  1.5  0.0  1.5
      1.5  0.0  0.0  0.0  1.5  0.0
   ]
   @test Oxygen.valencies(benzene) == [2, 2, 2, 2, 2, 2]
   @test Oxygen.atomtypes(benzene) == [
      1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
      1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
      1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
      1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
      1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
      1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
   ]
end

@testset "AtomOps" begin
   ethyl_acetate = smilestomol("CC(=O)OCC")
   @test map(Oxygen.valence_electrons, ethyl_acetate.atoms) == [4, 4, 6, 6, 4, 4]
   @test Oxygen.neighbors(ethyl_acetate, 2) == [OxygenAtom("C"), OxygenAtom("O"), OxygenAtom("O")]

   benzene = smilestomol("c1ccccc1")
   @test map(Oxygen.valence_electrons, benzene.atoms) == [4, 4, 4, 4, 4, 4]
   @test Oxygen.neighbors(benzene, 2) == [OxygenAtom("C"), OxygenAtom("C")]
end
