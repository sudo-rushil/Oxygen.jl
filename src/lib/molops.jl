#= Mol Operations
Functions to process raw OxygenMols and extract tangential low level information from such.

=#

function process!(mol::OxygenMol)::OxygenMol
    for (atom, bonds) in zip(mol.atoms, mol.adj)
        atom.degree = length(bonds)
    end
    mol
end


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


function getadjmatrix(mol::OxygenMol, ::Type{T})::Array{T,2} where T<:AbstractFloat
    adj = zeros(T, (length(mol.atoms),length(mol.atoms)))
    for (start, bonds) in enumerate(mol.adj)
        for (last, type) in bonds
            adj[start,last] = type
        end
    end
    adj
 end

 getadjmatrix(mol::OxygenMol) = getadjmatrix(mol::OxygenMol, Float32)
