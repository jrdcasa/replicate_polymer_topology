import os
import re
import numpy as np
from _collections import defaultdict


def insert_nonopenmm_ff_terms(fnamepdb, non_openmm_potentials_terms):

    base, _ = os.path.splitext(fnamepdb)
    base = os.path.split(base)[-1]
    dirlocal = "./"

    # Get all content of the file
    ext = ".top"
    # cJ fname = os.path.join(dirlocal, base + "_replicate" + ext)
    fname = os.path.join(dirlocal, base + ext)
    with open(fname, "r") as f:
        contents = f.readlines()

    # Locate moleculetype labels
    idx_moleculetype_labels = []
    idx_dihedral_insert = []
    idx_insert_return_atend = []
    idx_atoms_labels = []
    for idx in range(len(contents)):
        if contents[idx] == "[ moleculetype ]\n":
            idx_moleculetype_labels.append(idx)
        if contents[idx] == "[ atoms ]\n":
            idx_atoms_labels.append(idx)
    idx_moleculetype_labels.append(len(contents))

    # Need map global number to local numbering for each kind of molecules
    nkinds = len(non_openmm_potentials_terms['toxwaerd']["torsiononlyonekind"])
    iglobal_at = 0
    maps_globatidx_to_localidx = defaultdict()
    ikind_idx = 0
    kinds_labels = list(non_openmm_potentials_terms['toxwaerd']["torsiononlyonekind"].keys())
    if nkinds >= 1:
        for idx_atom in idx_atoms_labels:
            j = idx_atom + 1
            while True:
                if len(contents[j]) < 2:
                    ikind_idx += 1
                    try:
                        iglobal_at = non_openmm_potentials_terms["toxwaerd"]["atoms_first_molecule_each_kind_dict"]\
                                                                [kinds_labels[ikind_idx]][0]
                    except IndexError:
                        pass
                    break
                elif re.match("^;", contents[j]):
                    j += 1
                else:
                    iline = contents[j]
                    maps_globatidx_to_localidx[iglobal_at+1] = int(iline.split()[0])
                    j += 1
                    iglobal_at += 1

    # Locate point to insert dihedrals
    for idx in range(len(idx_moleculetype_labels)-1):
        istart = idx_moleculetype_labels[idx]
        iend = idx_moleculetype_labels[idx+1]
        # Find the position to insert the new dihedrals
        pos = -1
        # Insert after [ dihedrals ] labels, if exists
        j = istart
        for iline in contents[istart:iend]:
            if iline == "[ dihedrals ]\n":
                pos = j + 1
                idx_insert_return_atend.append(False)
                break
            else:
                j += 1
        # Insert after [ angles ] labels.
        if pos == -1:
            j = istart
            for iline in contents[istart:iend]:
                if iline == "[ angles ]\n":
                    pos = j
                    for jline in contents[pos:iend]:
                        if len(jline) < 2:
                            pos = j + 1
                            # contents.insert(pos, "[ dihedrals ]\n")
                            idx_dihedral_insert.append(pos)
                            idx_insert_return_atend.append(True)
                            break
                        else:
                            j += 1
                else:
                    j += 1
        else:
            idx_dihedral_insert.append(pos)

    # Produce the lines
    for key in non_openmm_potentials_terms:
        if key == "toxwaerd":
            toxward_lines = setup_onlyonemol_toxwaerd_lines(non_openmm_potentials_terms['toxwaerd'],
                                                            maps_globatidx_to_localidx)

    i = 0
    for ikey, ivalues in toxward_lines.items():
        if idx_insert_return_atend[i]:
            l = "\n[ dihedrals ]\n"
            ivalues = l + ivalues
            ivalues += "\n"
            contents.insert(idx_dihedral_insert[i], ivalues)
            for j in range(i+1, len(idx_dihedral_insert)):
                idx_dihedral_insert[j] = idx_dihedral_insert[j] + 1
        else:
            contents.insert(idx_dihedral_insert[i], ivalues)
            for j in range(i+1, len(idx_dihedral_insert)):
                idx_dihedral_insert[j] = idx_dihedral_insert[j] + 1
        i += 1

    with open(fname, "w") as f:
        contents = "".join(contents)
        f.write(contents)


# =============================================================================
def setup_all_toxwaerd_lines(toxwaerd_dict):

    ilines = ";multiple propers (Toxwaerd torsion)\n"

    for idih in toxwaerd_dict["torsions"]:
        atoms_dih = list(idih[0])
        type_dih_forward = idih[1]
        if type_dih_forward in toxwaerd_dict.keys():
            type_dih = type_dih_forward
        else:
            ll = idih[1].split("-")
            type_dih_reverse = ll[3] + "-" + ll[2] + "-" + ll[1] + "-" + ll[0]
            type_dih = type_dih_reverse
        for idx in range(0, 9):
            ilines += "{0:6d} {1:6d} {2:6d} {3:6d} 9 {4:5.1f} {5:7.4f} {6:1d} ;\n"\
                .format(atoms_dih[0]+1, atoms_dih[1]+1,
                        atoms_dih[2]+1, atoms_dih[3]+1,
                        float(toxwaerd_dict[type_dih]["phase"][idx])*180./np.pi,
                        float(toxwaerd_dict[type_dih]["coeff"][idx]),
                        int(toxwaerd_dict[type_dih]["periodicity"][idx]))

    return ilines


# =============================================================================
def setup_onlyonemol_toxwaerd_lines(toxwaerd_dict, maps_globatidx_to_localidx):

    ilines = defaultdict(str)
    nkinds = len(toxwaerd_dict["torsiononlyonekind"])

    for ikind, dihedrals in toxwaerd_dict["torsiononlyonekind"].items():
        ilines[ikind] = ";multiple propers (Toxwaerd torsion)\n"

        for idih in dihedrals:
            atoms_dih = list(idih[0])
            type_dih_forward = idih[1]
            if type_dih_forward in toxwaerd_dict.keys():
                type_dih = type_dih_forward
            else:
                ll = idih[1].split("-")
                type_dih_reverse = ll[3] + "-" + ll[2] + "-" + ll[1] + "-" + ll[0]
                type_dih = type_dih_reverse
            # If not the dihedral have already be in the topology file
            if type_dih in toxwaerd_dict.keys():
                for idx in range(0, 9):
                    i = atoms_dih[0]+1
                    j = atoms_dih[1]+1
                    k = atoms_dih[2]+1
                    l = atoms_dih[3]+1
                    iloc1 = maps_globatidx_to_localidx[i]
                    iloc2 = maps_globatidx_to_localidx[j]
                    iloc3 = maps_globatidx_to_localidx[k]
                    iloc4 = maps_globatidx_to_localidx[l]
                    ilines[ikind] += "{0:6d} {1:6d} {2:6d} {3:6d} 9 {4:5.1f} {5:7.4f} {6:1d} ;\n"\
                        .format(iloc1, iloc2,
                                iloc3, iloc4,
                                float(toxwaerd_dict[type_dih]["phase"][idx])*180./np.pi,
                                float(toxwaerd_dict[type_dih]["coeff"][idx]),
                                int(toxwaerd_dict[type_dih]["periodicity"][idx]))

    return ilines
