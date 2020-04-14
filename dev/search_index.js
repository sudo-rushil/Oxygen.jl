var documenterSearchIndex = {"docs":
[{"location":"todo.html#TODO-1","page":"TODO","title":"TODO","text":"","category":"section"},{"location":"todo.html#","page":"TODO","title":"TODO","text":"As Oxygen.jl is in active development, there are many things to be done. Any contributions would be greatly appreciated.","category":"page"},{"location":"todo.html#","page":"TODO","title":"TODO","text":"If you want to help, here are some ideas:","category":"page"},{"location":"todo.html#","page":"TODO","title":"TODO","text":"Help write tests or documentation.\nAdd functionality (see project structure below.)","category":"page"},{"location":"todo.html#Project-Structure-1","page":"TODO","title":"Project Structure","text":"","category":"section"},{"location":"todo.html#","page":"TODO","title":"TODO","text":"This list details the expected eventual structure for the Oxygen.jl package. Use this as a guide to the file heirarchy. Keep in mind not all these files exist yet.","category":"page"},{"location":"todo.html#","page":"TODO","title":"TODO","text":"Chem\nPeriodic Table - Basic elemental constants.\nElectrons - Electron configuration types.\nAtoms - Atom types.\nMolecules - Graph-based organic molecule types.\nReactions - Wrapper types for organic reactions.\nLib\nMolOps - Molecular processing and data extraction.\nDescriptors - Molecular descriptors functions (steric, surface area, etc.).\nFingerprints - Molecular fingerprints (ECFP, Morgan).\nForce Field - Functions for electronic force fields of molecules.\nCoordinates - 2D and 3D atomic coordinates.\nIO - Read and write OxygenMols from different formats.\nSMILES\nSMARTS\nInChi\nMolfile\nSDF\nCML\nQuery\nSubstructure - Find molecular substructural motifs in an OxygenMol.\nRing Finding - Find rings and aromatic rings in a molecule.\nDraw - Create structural formula from molecules.","category":"page"},{"location":"api/smiles.html#SMILES-Parser-1","page":"SMILES Parser","title":"SMILES Parser","text":"","category":"section"},{"location":"api/smiles.html#","page":"SMILES Parser","title":"SMILES Parser","text":"Functions to convert smiles strings into OxygenMol objects.","category":"page"},{"location":"api/smiles.html#","page":"SMILES Parser","title":"SMILES Parser","text":"SMILES parser is expected to change and evolve over time","category":"page"},{"location":"api/smiles.html#","page":"SMILES Parser","title":"SMILES Parser","text":"Modules = [Oxygen]\nPages   = [\"smiles.jl\"]\nOrder   = [:function]","category":"page"},{"location":"api/smiles.html#Oxygen.smilestomol-Tuple{String}","page":"SMILES Parser","title":"Oxygen.smilestomol","text":"smilestomol(smiles)\n\nReturns the OxygenMol for the molecule represented by smiles.\n\n\n\n\n\n","category":"method"},{"location":"api/molops.html#MolOps-1","page":"MolOps","title":"MolOps","text":"","category":"section"},{"location":"api/molops.html#","page":"MolOps","title":"MolOps","text":"Operations to process and extract information from OxygenMols.","category":"page"},{"location":"api/molops.html#","page":"MolOps","title":"MolOps","text":"Modules = [Oxygen]\nPages   = [\"molops.jl\"]\nOrder   = [:function]","category":"page"},{"location":"index.html#Oxygen.jl-Documentation-1","page":"Home","title":"Oxygen.jl Documentation","text":"","category":"section"},{"location":"index.html#","page":"Home","title":"Home","text":"Oxygen.jl is a package for cheminformatics. It is focused on providing efficient support for chemical structure analysis and deep learning over molecules, as well as general purpose chemical data manipulation in Julia.","category":"page"},{"location":"index.html#","page":"Home","title":"Home","text":"Keep in mind Oxygen.jl is under early stages of active development. Features are expected to change frequently.","category":"page"},{"location":"index.html#Introduction-1","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"index.html#","page":"Home","title":"Home","text":"This website serves as documentation for the Oxygen.jl package. The sections of the documentation are accessible on the left hand sidebar.","category":"page"},{"location":"index.html#","page":"Home","title":"Home","text":"If you want documentation for a particular function, please browse the API or use the search functionality.","category":"page"},{"location":"index.html#Installing-1","page":"Home","title":"Installing","text":"","category":"section"},{"location":"index.html#","page":"Home","title":"Home","text":"Currently, Oxygen.jl only offers limited functionality. If you want to install it, use Julia's package manager. Type Pkg.clone(\"https://github.com/sudo-rushil/Oxygen.jl.git\") at the Julia REPL. Then, it can be used like any other Julia package.","category":"page"},{"location":"index.html#","page":"Home","title":"Home","text":"using Oxygen","category":"page"},{"location":"index.html#Contributing-1","page":"Home","title":"Contributing","text":"","category":"section"},{"location":"index.html#","page":"Home","title":"Home","text":"Contributions would be very welcomed! For ideas of what would help, please check the TODO section of this documentation. Please see the repository on Github and feel free to open issues or submit pull requests!","category":"page"},{"location":"index.html#Author-1","page":"Home","title":"Author","text":"","category":"section"},{"location":"index.html#","page":"Home","title":"Home","text":"This package is being written by Rushil Mallarapu. If you find it helpful, shoot me an email!","category":"page"},{"location":"api/atoms.html#OxygenAtom-Type-1","page":"OxygenAtom Type","title":"OxygenAtom Type","text":"","category":"section"},{"location":"api/atoms.html#","page":"OxygenAtom Type","title":"OxygenAtom Type","text":"Mutable types for representing atoms and associated atom features.","category":"page"},{"location":"api/atoms.html#","page":"OxygenAtom Type","title":"OxygenAtom Type","text":"Modules = [Oxygen]\nPages   = [\"atom.jl\"]\nOrder = [:type]","category":"page"},{"location":"api/mols.html#OxygenMol-Type-1","page":"OxygenMol Type","title":"OxygenMol Type","text":"","category":"section"},{"location":"api/mols.html#","page":"OxygenMol Type","title":"OxygenMol Type","text":"Data types for organic molecules.","category":"page"},{"location":"api/mols.html#","page":"OxygenMol Type","title":"OxygenMol Type","text":"OxygenMol implements a mutable molecule type, using weighted adjacency lists to denote bonds. A bond order of 1.5 refers to an aromatic bond. Hydrogens are implicit.","category":"page"},{"location":"api/mols.html#","page":"OxygenMol Type","title":"OxygenMol Type","text":"Modules = [Oxygen]\nPages   = [\"molecules.jl\"]\nOrder = [:type]","category":"page"}]
}