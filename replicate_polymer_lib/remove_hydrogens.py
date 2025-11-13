from collections import defaultdict
import os
import warnings
from replicate_polymer_lib.check_connect_pdb import check_and_remove_ter_labels
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    import mdtraj as md
warnings.filterwarnings("ignore")


# =============================================================================
def get_tempfactor_occup_list(pdb):

    """
    Get the temp factor list from a pdb file

    Args:
        pdb:

    Returns:
    A list containing the temperature factors
    """

    tempfactor_list = []
    occupancy_list = []
    with open(pdb, 'r') as fpdb:
        lines = fpdb.readlines()
        for iline in lines:
            if iline.find("ATOM") != -1 or iline.find("HETATM") != -1:
                try:
                    tempfactor_list.append(float(iline[60:66]))
                except ValueError:
                    tempfactor_list.append(float(0.0))

                try:
                    occupancy_list.append(float(iline[54:60]))
                except ValueError:
                    occupancy_list.append(float(0.0))


    return tempfactor_list, occupancy_list


# =============================================================================
def remove_hydrogens(filename_pdb, bondlist=None, removeallh=False):
    """
    This function takes the name of a pdb file and removes the hydrogen atoms.
    This can be usefult oe be used in united group force fields.

    Args:
        filename_pdb(str): Name of the file containing the coordinates without atomic order.
        bondlist(list)
        removeallh (bool): If True remove all Hs, otherwise only remove Hs attached to carbon atoms (united atoms)

    Returns:
        Name of the new pdb file

    """

    # Read the original pdb in which the atoms are not ordered
    traj = md.load_pdb(filename_pdb)
    atoms_UA = list(traj.topology.select("element != H"))
    tempfactor_single, occup_single = get_tempfactor_occup_list(filename_pdb)

    # Create boinds if there are not present in the md.traj
    if bondlist is not None and traj.topology.n_bonds == 0:
        for ibond in bondlist:
            iat1_idx = ibond[0]
            iat2_idx = ibond[1]
            iat1 = traj.topology._atoms[iat1_idx]
            iat2 = traj.topology._atoms[iat2_idx]
            traj.topology.add_bond(iat1, iat2)

    # Add hydrogens explicitly to atoms different to C
    dict_graph = defaultdict(list)
    if not removeallh:
        for ibond in traj.topology.bonds:
            e1 = ibond[0].element.__str__()
            e2 = ibond[1].element.__str__()
            dict_graph[ibond[0].index].append(ibond[1].index)
            dict_graph[ibond[1].index].append(ibond[0].index)
            if e1 == "hydrogen" and e2 != "carbon":
                atoms_UA.append(ibond[0].index)
            if e2 == "hydrogen" and e1 != "carbon":
                atoms_UA.append(ibond[1].index)
        atoms_UA = sorted(atoms_UA)
    else:
        for ibond in traj.topology.bonds:
            e1 = ibond[0].element.__str__()
            e2 = ibond[1].element.__str__()
            dict_graph[ibond[0].index].append(ibond[1].index)
            dict_graph[ibond[1].index].append(ibond[0].index)

    tempfactor = []
    occfactor = []
    for iatom in atoms_UA:
        tempfactor.append(tempfactor_single[iatom])
        occfactor.append(occup_single[iatom])

    # Rename united atoms
    idx = 0
    atoms = [i for i in traj.topology.atoms]
    for iatom in atoms:
        nHs = 0
        nHeavy = 0
        if iatom.element.__str__() == "carbon":
            for jdx in dict_graph[idx]:
                if atoms[jdx].element.__str__() == "hydrogen":
                    nHs += 1
                elif atoms[jdx].element.__str__() != "carbon":
                    nHeavy += 1
            if nHs == 3:
                iatom.name = "_CH3"
            elif nHs == 2:
                iatom.name = "_CH2"
            elif nHs == 1:
                iatom.name = "_CH"
            elif nHs == 0:
                iatom.name = "C"
            elif nHs == 4:
                iatom.name = "_CH4"

        idx += 1

    new_trj = traj.atom_slice(atoms_UA)

    # Renumber the atoms in the new trajectory
    for ibond in new_trj.topology.bonds:
        a1 = ibond[0]
        a2 = ibond[1]
        a1.serial = a1.index + 1
        a2.serial = a2.index + 1

    # PDB
    baset, extt = os.path.splitext(os.path.basename(filename_pdb))
    newname_pdb = baset + "_noH" + extt
    new_trj.save_pdb(newname_pdb, bfactors=tempfactor)
    filenamepdb = check_and_remove_ter_labels(newname_pdb, occup=occfactor, logger=None)

    # GRO
    extt = ".gro"
    newname_gro = baset + "_noH" + extt
    new_trj.save_gro(newname_gro)

    # Remove TER labels


    return newname_pdb, newname_gro

