#= Atom Operations
Functions to extract information from OxygenAtoms.

=#

function valence_electrons(a::OxygenAtom)::Int
    sum(map(e -> e[end], a.orbitals.orbitals))
end


function neighbors(mol::OxygenMol, idx::Int)::Array{OxygenAtom,1}
    atom_neighbors = mol.adj[idx]
    [mol.atoms[i] for (i, _) in atom_neighbors]
end
