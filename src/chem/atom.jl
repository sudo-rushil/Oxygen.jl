#= Atom
Implements static and mutable types for representing atoms and associated atom features.

StaticAtom implements properties that do not change depending on atomic context, but OxygenAtom contains properties that do.

=#

import Base

abstract type Atomic end

struct StaticAtom <: Atomic
    number::Int
    symbol::String
    mass::Float64
    orbitals::Electrons
    electronegativity::Union{Float64,Missing}
end

mutable struct OxygenAtom <: Atomic
    number::Int
    symbol::String
    mass::Float64
    orbitals::Electrons
    electronegativity::Float64
    formalcharge::Int
    radicals::Int
    degree::Int
end

StaticAtom(symbol::String) = begin
    number = periodic_table[symbol]
    StaticAtom(number, symbol, periodic_table_mass[number], Electrons(number), periodic_table_electronegativity[number])
end

StaticAtom(number::Int) = begin
    symbol = periodic_table_reverse[number]
    StaticAtom(number, symbol, periodic_table_mass[number], Electrons(number), periodic_table_electronegativity[number])
end

OxygenAtom(number, symbol, mass, orbitals, electronegativity) = OxygenAtom(number, symbol, mass, orbitals, electronegativity, 0, 0, 0)

OxygenAtom(symbol::String) = begin
    number = periodic_table[symbol]
    OxygenAtom(number, symbol, periodic_table_mass[number], Electrons(number), periodic_table_electronegativity[number])
end

OxygenAtom(number::Int) = begin
    symbol = periodic_table_reverse[number]
    OxygenAtom(number, symbol, periodic_table_mass[number], Electrons(number), periodic_table_electronegativity[number])
end

function atom_equal(a::Atomic, b::Atomic)::Bool
    a.number == b.number && a.symbol == b.symbol && a.mass == b.mass
end

Base.show(io::IO, a::Atomic) = print(io, "\"", a.symbol, "\" atom (number ", a.number, ")")

# Base.show(io::IO, ::MIME"text/plain", a::Atomic) =
#            print(io, "$(typeof(a)):\n   ", a)
