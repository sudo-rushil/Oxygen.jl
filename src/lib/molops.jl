#= Mol Operations
Functions to process raw OxygenMols and extract tangential low level information from such.

=#


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


function atomtypes(::Type{T}, mol::OxygenMol)::Array{T,2} where {T<:AbstractFloat}
    validatoms = String["C", "N", "O", "S", "F", "P", "Cl", "Mg", "Na", "Br", "Ca", "Cu", "Pd", "I", "Al"]
    atoms = map(a -> _onehot(validatoms, a.symbol), mol.atoms)
    transpose(hcat(atoms...))
end

atomtypes(mol::OxygenMol) = atomtypes(Float32, mol::OxygenMol)


function valencies(mol::OxygenMol)::Array{Int,2}
    valcount = [length(bonds) for bonds in mol.adj]
    valences = map(a -> _onehot(collect(1:6), a), valcount)
    transpose(hcat(valences...))
end


function hybridization(mol::OxygenMol)::Array{Int,2}
    #TODO: Retrieve hybridization states instead of bonding electrons.
    electrons(bonds) = sum([e for (_, e) in bonds])
    hybcount = [electrons(bonds) for bonds in mol.adj]
    hybrid = map(a -> _onehot(collect(1:.5:4), a), hybcount)
    transpose(hcat(hybrid...))
end


function molfeatures(::Type{T}, mol::OxygenMol)::Array{T, 2} where {T <: AbstractFloat}
    atoms = atomtypes(T, mol)
    valence = convert(Array{T}, valencies(mol))
    hybrid = convert(Array{T}, hybridization(mol))
    hcat(atoms, valence, hybrid)
end

molfeatures(mol::OxygenMol) = molfeatures(Float32, mol::OxygenMol)
