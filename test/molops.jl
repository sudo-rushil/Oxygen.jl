using Oxygen

println("Testing Mol Operations")

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
end
