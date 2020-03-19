#= Molecules
Data types for organic molecules.

OxygenMol implements a mutable molecule type, using weighted adjacency lists to denote bonds.
A bond order of 1.5 refers to an aromatic bond.
Hydrogens are implicit

TODO: Implement a static molecule type and improve equality testing
=#

abstract type Molecular end

mutable struct OxygenMol <: Molecular
    atoms::Array{OxygenAtom,1}
    adj::Array{Array{Tuple{Int,Int},1},1}
end

OxygenMol(atoms::Array{Int,1}) = OxygenMol(map(number -> OxygenAtom(number), atoms), [])

OxygenMol(atoms::Array{String,1}) = OxygenMol(map(symbol -> OxygenAtom(symbol), atoms), [])

function mol_equal(a::Molecular, b::Molecular)::Bool
    a.adj == b.adj
end
