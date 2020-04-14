# TODO

As Oxygen.jl is in active development, there are many things to be done. Any contributions would be greatly appreciated.

If you want to help, here are some ideas:
- Help write tests or documentation.
- Add functionality (see project structure below.)

# Project Structure

This list details the expected eventual structure for the Oxygen.jl package. Use this as a guide to the file heirarchy. Keep in mind not all these files exist yet.

- Chem
  - Periodic Table - Basic elemental constants.
  - Electrons - Electron configuration types.
  - Atoms - Atom types.
  - Molecules - Graph-based organic molecule types.
  - Reactions - Wrapper types for organic reactions.
- Lib
  - MolOps - Molecular processing and data extraction.
  - Descriptors - Molecular descriptors functions (steric, surface area, etc.).
  - Fingerprints - Molecular fingerprints (ECFP, Morgan).
  - Force Field - Functions for electronic force fields of molecules.
  - Coordinates - 2D and 3D atomic coordinates.
- IO - Read and write OxygenMols from different formats.
  - SMILES
  - SMARTS
  - InChi
  - Molfile
  - SDF
  - CML
- Query
  - Substructure - Find molecular substructural motifs in an OxygenMol.
  - Ring Finding - Find rings and aromatic rings in a molecule.
- Draw - Create structural formula from molecules.
