import os
import numpy as np


def insert_nonopenmm_ff_terms(fnamepdb, non_openmm_potentials_terms):

    base, _ = os.path.splitext(fnamepdb)
    base = os.path.split(base)[-1]
    dirlocal = "./"

    # Get all content of the file
    ext = ".top"
    fname = os.path.join(dirlocal, base + "_replicate" + ext)
    with open(fname, "r") as f:
        contents = f.readlines()

    # Find the position to insert the new dihedrals
    idx_to_insert = contents.index("[ dihedrals ]\n")
    for iline in contents[idx_to_insert:]:
        if len(iline) < 2:
            break
        else:
            idx_to_insert += 1

    # Produce the lines
    for key in non_openmm_potentials_terms:
        if key == "toxwaerd":
            toxward_lines = setup_all_toxwaerd_lines(non_openmm_potentials_terms['toxwaerd'])

    contents.insert(idx_to_insert, toxward_lines)

    with open(fname, "w") as f:
        contents = "".join(contents)
        f.write(contents)


# =============================================================================
def setup_all_toxwaerd_lines(toxwaerd_dict):

    ilines = ";multiple propers (Toxwaerd torsion)\n"

    for idih in toxwaerd_dict["torsions"]:
        atoms_dih = list(idih[0])
        type_dih = idih[1]
        for idx in range(0, 9):
            ilines += "{0:6d} {1:6d} {2:6d} {3:6d} 9 {4:5.1f} {5:7.4f} {6:1d} ;\n"\
                .format(atoms_dih[0]+1, atoms_dih[1]+1,
                        atoms_dih[2]+1, atoms_dih[3]+1,
                        float(toxwaerd_dict[type_dih]["phase"][idx])*180./np.pi,
                        float(toxwaerd_dict[type_dih]["coeff"][idx]),
                        int(toxwaerd_dict[type_dih]["periodicity"][idx]))

    return ilines