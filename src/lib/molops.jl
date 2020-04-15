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
