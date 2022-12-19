# Searching the Cambridge Structural Database for Relevant Crystal Structures

Purpose - showcase a method of searching for reported palladium complexes in the Cambridge Structural Database and then comparing to free ligand conformations.

```mermaid
graph TD;
    A(Substructure search of phosphine ligand SMILES using CSD API) -->|get_pr3_3d_coords.py| B(Phosphine-ligated palladium complex crystal structures);
    B-->|rdkit_check_pr3.py| C(3D coordinates of phosphine ligands of interest, confirm ligand match);
    C-->D(Consider conformation changes between Pd-ligated, phosphines, and phosphine oxides);
    D-->E(Xtal PR3);
    D-->F(Molecular mechanics PR3);
    D-->G(Molecular mechanics OPR3);
    E<-->F
    F<-->G
```
