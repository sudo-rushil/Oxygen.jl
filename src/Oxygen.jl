module Oxygen

export StaticAtom, OxygenAtom, OxygenMol, smilestomol, molfeatures

include("chem/periodic.jl")
include("chem/electrons.jl")
include("chem/atom.jl")
include("chem/molecules.jl")

include("lib/atomops.jl")
include("lib/molops.jl")

include("io/smiles.jl")

end
