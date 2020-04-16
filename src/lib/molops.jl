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


function valencies(mol::OxygenMol)::Array{Int,1}
    [length(bonds) for bonds in mol.adj]
end


function _onehot(states::Array{T,1}, state::T)::Array{Float32,1} where {T}
    [Float32(state == s) for s in states]
end


function atomtypes(::Type{T}, mol::OxygenMol)::Array{T, 2} where {T<:AbstractFloat}
    validatoms = String["C", "N", "O", "S", "F", "P", "Cl", "Mg", "Na", "Br", "Ca", "Cu", "Pd", "I", "Al"]
    atoms = map(a -> _onehot(validatoms, a.symbol), mol.atoms)
    transpose(hcat(atoms...))
end

atomtypes(mol::OxygenMol) = atomtypes(Float32, mol::OxygenMol)
