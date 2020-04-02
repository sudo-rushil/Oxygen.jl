using Oxygen

println("Testing SMILES Parser")

function sameatom(atoms::Tuple{OxygenAtom,OxygenAtom})::Bool
    (a,b) = atoms
    a.number == b.number
end

function sameatoms(a::Array{OxygenAtom,1}, b::Array{OxygenAtom,1})::Bool
    all(map(sameatom, zip(a,b)))
end

@testset "SmilesToMol" begin
    ethyl_acetate = "CC(=O)OCC"
    ethyl_acetate_mol = smilestomol(ethyl_acetate)
    @test typeof(ethyl_acetate_mol) == OxygenMol
    @test sameatoms(ethyl_acetate_mol.atoms, [OxygenAtom("C"), OxygenAtom("C"), OxygenAtom("O"), OxygenAtom("O"), OxygenAtom("C"), OxygenAtom("C")])
    @test ethyl_acetate_mol.adj ==
          [[(2, 1.0)], [(1, 1.0), (3, 2.0), (4, 1.0)], [(2, 2.0)], [(2, 1.0), (5, 1.0)], [(4, 1.0), (6, 1.0)], [(5, 1.0)]]

    benzene = "c1ccccc1"
    benzene_mol = smilestomol(benzene)
    @test typeof(benzene_mol) == OxygenMol
    @test sameatoms(benzene_mol.atoms, [OxygenAtom("C"), OxygenAtom("C"), OxygenAtom("C"), OxygenAtom("C"), OxygenAtom("C"), OxygenAtom("C")])
    @test benzene_mol.adj == [
        [(2, 1.5), (6, 1.5)],
        [(1, 1.5), (3, 1.5)],
        [(2, 1.5), (4, 1.5)],
        [(3, 1.5), (5, 1.5)],
        [(4, 1.5), (6, 1.5)],
        [(5, 1.5), (1, 1.5)],
    ]

    cipro = "C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O"
    cipro_mol = smilestomol("C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O")
    @test typeof(cipro_mol) == OxygenMol
    @test sameatoms(cipro_mol.atoms[7:11], [OxygenAtom("C"), OxygenAtom("O"), OxygenAtom("C"), OxygenAtom("C"), OxygenAtom("C")])
    @test cipro_mol.adj[7:11] == [
        [(6, 1.0), (8, 2.0), (9, 1.0)],
        [(7, 2.0)],
        [(7, 1.0), (10, 2.0), (14, 1.0)],
        [(9, 2.0), (11, 1.0)],
        [(10, 1.0), (12, 2.0), (21, 1.0)],
    ]
end
