#= SMILES IO
Functions to convert smiles strings into OxygenMol objects.

SMILES parser is expected to change and evolve over time

Grammar specification based on https://biocyc.org/help.html?object=smiles

=#

function _split_atoms(smiles::String)::OxygenMol
    atom_match = r"b|c|n|o|p|s|B|C|N|O|P|S|F|Cl|Br|I"
    OxygenMol(map(x -> uppercasefirst(String(x.match)), collect(eachmatch(atom_match, smiles))))
end

mutable struct _ParseNode
    first::Int
    last::Int
    children::Array{_ParseNode,1}
end

_ParseNode(f::Int, l::Int) = _ParseNode(f, l, [])

function _addchild!(node::_ParseNode, child::_ParseNode)
    push!(node.children, child)
    node
end

function _split_smiles(smiles::String)::Array{Any,1}
    splitsmiles = map(x -> String(x.match), collect(eachmatch(r"b|c|n|o|p|s|B|C|N|O|P|S|F|Cl|Br|I|\(|\)", smiles)))
    counter, current = 1, 1
    out = []

    for char in splitsmiles
        if (char == "(" || char == ")")
            if current != counter
                push!(out, _ParseNode(current, counter - 1))
                current = counter
            end
            push!(out, char)
        else
            counter += 1
        end
    end

    push!(out, _ParseNode(current, counter - 1))
end

function _parse_brackets(raw_branches::Array{Any,1})::_ParseNode
    last, lastchar = raw_branches[1], raw_branches[1]
    openbrackets = []
    closebrackets = []

    for branch in raw_branches
        if !(branch in ["(", ")"])
            if isempty(closebrackets) && !isempty(openbrackets)
                _addchild!(openbrackets[end], branch)
            elseif !isempty(closebrackets)
                _addchild!(closebrackets[end], branch)
            end
        end

        if branch == "("
            if lastchar != ")"
                closebrackets = []
                push!(openbrackets, last)
            end
        elseif branch == ")"
            if !isempty(openbrackets)
                push!(closebrackets, pop!(openbrackets))
            end
        else
            last = branch
        end

        lastchar = branch
    end
    raw_branches[1]
end

function _addbond!(mol::OxygenMol, first::Int, last::Int, type::Float64 = 1.0)::OxygenMol
    push!(mol.adj[first], (last, type))
    push!(mol.adj[last], (first, type))
    mol
end

function _split_bonds(smiles::String)::Array{String,1}
    map(x -> String(x.match), collect(eachmatch(r"(=|#)?(b|c|n|o|p|s|B|C|N|O|P|S|F|Cl|Br|I)", smiles)))
end

function _bondvalue(bond_symbol::String)::Float64
    if bond_symbol == "="
        2
    elseif bond_symbol == "#"
        3
    elseif islowercase(bond_symbol[1])
        1.5
    else
        1
    end
end

function _bond_type(index::Int, bonds::Array{String,1})::Float64
    bond = bonds[index][1:1]
    _bondvalue(bond)
end

function _traverse!(mol::OxygenMol, bonds::Array{String,1}, root::_ParseNode)
    _nodetobonds!(mol, bonds, root)

    for node in root.children
        _addbond!(mol, root.last, node.first, _bond_type(node.first, bonds))
        _traverse!(mol, bonds, node)
    end
end

function _nodetobonds!(mol::OxygenMol, bonds::Array{String,1}, node::_ParseNode)::OxygenMol
    for i = node.first:node.last-1
        _addbond!(mol, i, i + 1, _bond_type(i + 1, bonds))
    end
    mol
end

function _addcycles!(mol::OxygenMol, bonds::Array{String,1}, smiles::String)
    ring_atoms = map(x -> String(x.match), collect(eachmatch(r"(C|c|N|n|O|F|S|Cl|Br|P)(=|#)?\d?", smiles)))
    rings = Dict{Char,Int}()

    for i = 1:length(ring_atoms)
        index = ring_atoms[i][end]
        if isdigit(index)
            if !(index in keys(rings))
                rings[index] = i
            else
                bond_symbol = ring_atoms[1][end-1:end-1]
                bond = _bondvalue(bond_symbol)
                _addbond!(mol, rings[index], i, bond)
            end
        end
    end

    mol
end

"""
    smilestomol(smiles)

Returns the OxygenMol for the molecule represented by `smiles`.
"""
function smilestomol(smiles::String)::OxygenMol
    # Temporary workaround for atom-labelling or explicit hydrogens.
    cleansmiles = replace(smiles, r"\[(b|c|n|o|p|s|B|C|N|O|P|S|F|Cl|Br|I).+?\]" => s"\1")

    mol = _split_atoms(cleansmiles)
    branches = _split_smiles(cleansmiles)
    bonds = _split_bonds(cleansmiles)
    root = _parse_brackets(branches)

    _traverse!(mol, bonds, root)
    _addcycles!(mol, bonds, cleansmiles)
    process!(mol)

    mol
end
