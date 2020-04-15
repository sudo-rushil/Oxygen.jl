#= Atom
Implements static and mutable types for representing atoms and associated atom features.

Depreciated: StaticAtom implements properties that do not change depending on atomic context, but OxygenAtom contains properties that do.
Current: OxygenAtom is a static type, and other variable properties are always calculated with respect to their environment.

=#

abstract type Atomic end

struct OxygenAtom <: Atomic
    number::Int
    symbol::String
    mass::Float64
    orbitals::Electrons
    electronegativity::Float64
end

OxygenAtom(symbol::String) = begin
    number = periodic_table[symbol]
    OxygenAtom(number, symbol, periodic_table_mass[number], Electrons(number), periodic_table_electronegativity[number])
end

OxygenAtom(number::Int) = begin
    symbol = periodic_table_reverse[number]
    OxygenAtom(number, symbol, periodic_table_mass[number], Electrons(number), periodic_table_electronegativity[number])
end


Base.:(==)(a::Atomic, b::Atomic) = isequal(a.number, b.number) && isequal(a.symbol, b.symbol) && isequal(a.mass, b.mass) && isequal(a.orbitals, b.orbitals)

Base.show(io::IO, a::Atomic) = print(io, "\"", a.symbol, "\" atom (number ", a.number, ")")

# Base.show(io::IO, ::MIME"text/plain", a::Atomic) =
#            print(io, "$(typeof(a)):\n   ", a)
