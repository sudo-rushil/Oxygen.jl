#= Molecules
Data types for organic molecules.

OxygenMol implements a mutable molecule type, using weighted adjacency lists to denote bonds.
A bond order of 1.5 refers to an aromatic bond.
Hydrogens are implicit

TODO: Implement a static molecule type
=#

import AutoHashEquals: @auto_hash_equals

abstract type Molecular end

@auto_hash_equals mutable struct OxygenMol <: Molecular
    atoms::Array{OxygenAtom,1}
    adj::Array{Array{Tuple{Int,Float64},1},1}
end

"""
    OxygenMol

Mutable molecule type.

Represents molecules using an array of atoms and a weighted adjacency list to represent bonds.
"""
OxygenMol

OxygenMol(atoms::Array{Int,1}) = OxygenMol(map(number -> OxygenAtom(number), atoms), [[] for _ in 1:length(atoms)])

OxygenMol(atoms::Array{String,1}) = OxygenMol(map(symbol -> OxygenAtom(symbol), atoms), [[] for _ in 1:length(atoms)])
