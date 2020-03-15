#= Electronic Structure
Implements types for dealing with electronic configurations of atoms.

TODO: See if tuple implementation of orbitals needs to be improved on.

=#

import Base

struct Electrons
    root::Union{String, Nothing}
    orbitals::Tuple
end

function parse_config(config::String)::Electrons
    config_pattern = r"(?:\[(\w+)\])?((?:\w|\s?)+)"
    captures = match(config_pattern, config).captures
    orbitals = sort(split(captures[2]), by = o -> o[1])

    Electrons(
        captures[1],
        tuple(
            map(o -> tuple(o[1:2], parse(Int, o[3:end])), orbitals)...
        )
    )
end

parse_config("[Ne]3s2 3p4")

Electrons(number::Int) = begin
    config = periodic_table_electrons[number]
    parse_config(config)
end

Electrons(symbol::String) =  begin
    number = periodic_table[symbol]
    Electrons(number)
end

Base.show(io::IO, e::Electrons) = print(io, "[", e.root, "] ", e.orbitals)

# Base.show(io::IO, ::MIME"text/plain", e::Electrons) =
#            print(io, "$(typeof(a)):\n   ", e)
