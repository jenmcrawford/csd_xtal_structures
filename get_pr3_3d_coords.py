from ccdc import io
import ccdc.search

def break_pd_bonds(mol2_file):
    """
    Function that breaks all bonds to palladium
    arg1: mol2 file containing Pd phosphine complexes to search
    Returns: list of CSD identifiers and nested list of molecules and fragments
    """
    mols = io.MoleculeReader(mol2_file)
    names = [m.identifier for m in mols]
    new_mols = []
    for m in mols:
        m.remove_bonds(b for b in m.bonds if any(a.atomic_number == 46 for a in b.atoms))
        new_mols.append(m.components)
    return names, new_mols

def is_pr3(fragments):
    """
    Function that identifies if a particular component contains phosphorus
    arg1: fragments of each molecule
    Returns: nested list of mols and fragments that contain phosphorus, needs to be further tested in RDKit (see script X)
    """

    p_mol = [[i for i in f for a in i.atoms if a.atomic_number == 15] for f in fragments]
    return p_mol


with open('phosphines_tosearch.txt') as f:
    full_list = f.readlines()
    gsk_id = [i.split(',',1)[0] for i in full_list]
    smi = [i.split(',',1)[1].strip() for i in full_list]

for x in range(len(gsk_id)):

    csd_id, fragments = break_pd_bonds(gsk_id[x] + '_pd_hits.mol2')
    p_mol = is_pr3(fragments)

    if len(p_mol) == len(csd_id):
        print("For " + gsk_id[x] + " the number of molecules match")

    for i in range(len(p_mol)):
        p_num = list(range(len(p_mol[i])))
        for j in p_num:
            with ccdc.io.MoleculeWriter(gsk_id[x] + '_' + csd_id[i] + '_' + str(j) + '.mol2') as writer:
                writer.write(p_mol[i][j])


