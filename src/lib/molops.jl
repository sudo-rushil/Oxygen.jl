#= Mol Operations
Functions to process raw OxygenMols and extract tangential low level information from such.

=#

"""
    getbonds(mol)

Retrieves list of bonds in a molecule.
Bond list formatted as `(starting atom, ending atom, bond type)`.
All bonds are unique.

# Examples
```
julia> getbonds(smilestomol("c1ccccc1"))
6-element Array{Tuple,1}:
 (1, 2, 1.5)
 (1, 6, 1.5)
 (2, 3, 1.5)
 (3, 4, 1.5)
 (4, 5, 1.5)
 (5, 6, 1.5)
```
"""
function getbonds(mol::OxygenMol)::Array{Tuple,1}
    bondlist = []
    for (start, bonds) in enumerate(mol.adj)
        for (last, type) in bonds
            if last > start
                push!(bondlist, (start, last, type))
            end
        end
    end
    bondlist
end

"""
    getadjmatrix(smiles)

Generates adjacency matrix of molecule, weighted by bond type

# Examples
```
julia> getadjmatrix(smilestomol("c1ccccc1"))
6×6 Array{Float32,2}:
 0.0  1.5  0.0  0.0  0.0  1.5
 1.5  0.0  1.5  0.0  0.0  0.0
 0.0  1.5  0.0  1.5  0.0  0.0
 0.0  0.0  1.5  0.0  1.5  0.0
 0.0  0.0  0.0  1.5  0.0  1.5
 1.5  0.0  0.0  0.0  1.5  0.0
```
"""
function getadjmatrix(::Type{T}, mol::OxygenMol)::Array{T,2} where {T<:AbstractFloat}
    adj = zeros(T, (length(mol.atoms), length(mol.atoms)))
    for (start, bonds) in enumerate(mol.adj)
        for (last, type) in bonds
            adj[start, last] = type
        end
    end
    adj
end

getadjmatrix(mol::OxygenMol) = getadjmatrix(Float32, mol::OxygenMol)


function _onehot(states::Array{T,1}, state::T)::Array{Float32,1} where {T}
    [Float32(state == s) for s in states]
end


"""
    atomtypes(smiles)

Creates one-hot mapping of atoms in a molecule to common atoms.
Current valid atoms array is `["C", "N", "O", "S", "F", "P", "Cl", "Mg", "Na", "Br", "Ca", "Cu", "Pd", "I", "Al"]`.

# Examples
```
julia> atomtypes(smilestomol("c1ccccc1"))
6×15 Array{Float32,2}:
 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
"""
function atomtypes(::Type{T}, mol::OxygenMol)::Array{T,2} where {T<:AbstractFloat}
    validatoms = String["C", "N", "O", "S", "F", "P", "Cl", "Mg", "Na", "Br", "Ca", "Cu", "Pd", "I", "Al"]
    atoms = map(a -> _onehot(validatoms, a.symbol), mol.atoms)
    transpose(hcat(atoms...))
end

atomtypes(mol::OxygenMol) = atomtypes(Float32, mol::OxygenMol)

"""
    valencies(smiles)

Creates one-hot mapping of atom valencies in a molecule.
Current maximum valid valencies (number of connected non-hydrogen atoms) is 6.

# Examples
```
julia> valencies(smilestomol("c1ccccc1"))
6×6 Array{Int64,2}:
 0  1  0  0  0  0
 0  1  0  0  0  0
 0  1  0  0  0  0
 0  1  0  0  0  0
 0  1  0  0  0  0
 0  1  0  0  0  0
```
"""
function valencies(mol::OxygenMol)::Array{Int,2}
    valcount = [length(bonds) for bonds in mol.adj]
    valences = map(a -> _onehot(collect(1:6), a), valcount)
    transpose(hcat(valences...))
end

"""
    hybridization(smiles)

Creates one-hot mapping of atom hybridization states in a molecule.
Currently creates onehot mapping of bonding electron count. Hybridization states will be added soon.

# Examples
```
julia> hybridization(smilestomol("c1ccccc1"))
6×7 Array{Int64,2}:
 0  0  0  0  1  0  0
 0  0  0  0  1  0  0
 0  0  0  0  1  0  0
 0  0  0  0  1  0  0
 0  0  0  0  1  0  0
 0  0  0  0  1  0  0
"""
function hybridization(mol::OxygenMol)::Array{Int,2}
    #TODO: Retrieve hybridization states instead of bonding electrons.
    electrons(bonds) = sum([e for (_, e) in bonds])
    hybcount = [electrons(bonds) for bonds in mol.adj]
    hybrid = map(a -> _onehot(collect(1:.5:4), a), hybcount)
    transpose(hcat(hybrid...))
end

"""
    molfeatures(smiles)

Creates atom-wise features matrix for a molecule based on atomtypes, valence, and hybridization.
Note: `molfeatures` and `getadjmatrix` do **not** return padded arrays.

# Examples
```
julia> molfeatures(smilestomol("c1ccccc1"))
6×28 Array{Float32,2}:
 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0
 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0
 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0
 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0
 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0
 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0
```
"""
function molfeatures(::Type{T}, mol::OxygenMol)::Array{T, 2} where {T <: AbstractFloat}
    atoms = atomtypes(T, mol)
    valence = convert(Array{T}, valencies(mol))
    hybrid = convert(Array{T}, hybridization(mol))
    hcat(atoms, valence, hybrid)
end

molfeatures(mol::OxygenMol) = molfeatures(Float32, mol::OxygenMol)
