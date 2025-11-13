from collections import defaultdict
import warnings
import datetime
import copy
import os
import re
import numpy as np
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    import mdtraj as md
from replicate_polymer_lib.remove_hydrogens import get_tempfactor_occup_list
from replicate_polymer_lib.check_connect_pdb import check_and_remove_ter_labels
import MDAnalysis as mdanalysis
warnings.filterwarnings("ignore")


# =============================================================================
def replicate_pdb(single_mol_pdb, images, boxlength=None, boxangles=None, save_single=False, logger_log=None):

    """

    Args:
        single_mol_pdb:
        images:
        boxlength:
        boxangles:
        save_single:
        logger_log:

    Returns:

    """

    # Convert pdb to mdtrj. Take some info from the pdb ==========
    # nch_single and trj_single.n_chains indicate the kind of molecules
    # mols_single is the total number of molecules
    trj_single = md.load_pdb(single_mol_pdb)
    topo_single = trj_single.topology
    try:
        mols_single = topo_single.find_molecules()
    except ValueError:
        _guessing_bonds(single_mol_pdb, topo_single, logger=logger_log)
        mols_single = topo_single.find_molecules()
    nmols_single = len(mols_single)
    nch_single = trj_single.n_chains
    nres_single = trj_single.n_residues
    nat_single = trj_single.n_atoms
    nframes_single = trj_single.n_frames
    nbonds_single = topo_single.n_bonds
    atom_to_chainid_single = _assign_chain_to_atom(single_mol_pdb)

    # MDTraj has several standardBonds for protein residues this can cause problems for non-protein systems in which
    # any residue has the same name that those standard residue names.
    standard_residues = topo_single._standardBonds.keys()
    for ires in topo_single.residues:
        if ires.name in standard_residues:
            mm1 = "\t\tERROR: The residue {} is a standard residue.\n".format(ires.name)
            mm2 = "\t\tThis can be a source of problems for non-protein systems\n"
            mm3 = "\t\tConsiderer to change the residue name before to use this program\n"
            mm0 = "\t\t"+len(mm3)*"-"+"\n"
            print(mm0+mm1+mm2+mm3+mm0) if logger_log is None else logger_log.info(mm0+mm1+mm2+mm3+mm0)
            exit()

    # Set of residues names in the PDB ==========
    name_residues = list(set([i.name for i in trj_single.topology.residues]))
    dict_nresidues = defaultdict(int)
    list_nresidues = []
    for ires in trj_single.topology.residues:
        iname = ires.name
        dict_nresidues[iname] += 1
        list_nresidues.append(iname)
    # Set of molecules names in the PDB ==========
    dict_nmolecules = defaultdict(list)
    for imol in mols_single:
        # Just check the first atom in each molecule. It is supposed that the rest belong to the same chainID(molecule)
        idx = list(imol)[0].index
        name = '%s%s' % ('system', atom_to_chainid_single[idx])
        try:
            dict_nmolecules[name][0] += 1
        except IndexError:
            dict_nmolecules[name].append(1)

    # NUmber of frames ==========
    if nframes_single > 1:
        mm = "\tNumber of frames in the original trajectory is greater than 1.\n"
        mm += "\t\tOnly the last frame will be used to create the replica."
        print(mm) if logger_log is None else logger_log.warning(mm)

    # Print info to log
    nows = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    mm1 = "\n\t\tReplicate a molecule or polymer chain ({})\n".format(nows)
    mm = "\t\t" + len(mm1) * "*" + "\n"
    mm += "\t\tInput file                  : {}\n".format(single_mol_pdb)
    mm += "\t\tNumber of molecules in cell : {}\n".format(nmols_single)
    mm += "\t\tNumber of chains in cell    : {}\n".format(nch_single)
    mm += "\t\tNumber of bonds in cell     : {}\n".format(nbonds_single)
    mm += "\t\tNumber of residues in cell  : {}\n".format(nres_single)
    mm += "\t\tNumber of atoms in cell     : {}\n".format(nat_single)

    # Box information --> If box info is not available, a guess of the lenghts is made.
    # else if CRYST info is available, the boxlength is overrided if the box_length option is present.
    if trj_single.unitcell_lengths is None and boxlength is None:
        xmin = min(trj_single.xyz[0, :, 0])
        ymin = min(trj_single.xyz[0, :, 1])
        zmin = min(trj_single.xyz[0, :, 2])
        for iatom in range(nat_single):
            trj_single.xyz[0, iatom, 0] -= xmin
            trj_single.xyz[0, iatom, 1] -= ymin
            trj_single.xyz[0, iatom, 2] -= zmin
        lx = max(trj_single.xyz[0, :, 0])
        ly = max(trj_single.xyz[0, :, 1])
        lz = max(trj_single.xyz[0, :, 2])
        boxlength = np.array([lx, ly, lz])
        trj_single.unitcell_lengths = boxlength
    elif boxlength is not None:
        trj_single.unitcell_lengths = boxlength
    else:
        boxlength = trj_single.unitcell_lengths[0, :]

    if trj_single.unitcell_angles is None and boxangles is None:
        boxangles = np.array([90, 90, 90])
        trj_single.unitcell_angles = boxangles
    elif boxangles is not None:
        trj_single.unitcell_angles = boxangles
    else:
        boxangles = trj_single.unitcell_angles[0, :]

    mm += "\t\tSimulation box dimension    : {0:.4f} {1:.4f} {2:.4f} nm\n".format(float(boxlength[0]),
                                                                                  float(boxlength[1]),
                                                                                  float(boxlength[2]))
    mm += "\t\tSimulation box angles       : {0:.2f} {1:.2f} {2:.2f} degrees\n".format(float(boxangles[0]),
                                                                                       float(boxangles[1]),
                                                                                       float(boxangles[2]))
    mm += "\t\tReplicate the molecule {}x{}x{} times\n".format(images[0], images[1], images[2])
    mm += "\t\t" + len(mm1) * "*" + "\n"
    print(mm1 + mm) if logger_log is None else logger_log.info(mm1 + mm)

    # Calculate the total number of atoms and molecules ==========
    nxx = images[0]
    nyy = images[1]
    nzz = images[2]
    nmols = nxx * nyy * nzz
    natoms = nmols * nat_single
    new_xyz = np.zeros([1, natoms, 3])

    # Initializate tempFactor list to take into account the backbone and branch atoms
    temp_factor_list_single, occup_list_single = get_tempfactor_occup_list(single_mol_pdb)
    tempFactor = np.zeros([natoms])
    occupFactor = np.zeros([natoms])

    # Copy the original cell ==========
    new_topo = copy.deepcopy(topo_single)
    mols_replicate = copy.deepcopy(mols_single)

    # If the original topology has not named the residues, give a name to each residue starting a 000 ==========
    ires = 0
    for _ in new_topo.chains:
        for residue in new_topo.residues:
            if residue.name == '':
                residue.name = "{0:03d}".format(ires)
                ires += 1
    new_topo._bonds = []
    for bond in topo_single.bonds:
        a1, a2 = bond
        a1n = new_topo.atom(a1.index)
        a2n = new_topo.atom(a2.index)
        new_topo.add_bond(a1n, a2n, bond.type, bond.order)

    # Coordinates of the first chain
    new_xyz[0, 0:nat_single, :] = trj_single.xyz[-1, :, :]
    tempFactor[0:nat_single] = temp_factor_list_single
    occupFactor[0:nat_single] = occup_list_single

    # Add the cells in the x-direction, y-direction and z-direction ==========
    ind_atom = nat_single
    resSeq = nres_single
    for inz in range(0, nzz):
        for iny in range(0, nyy):
            for inx in range(0, nxx):
                # First chain is already placed
                if inx == 0 and iny == 0 and inz == 0:
                    nreplicas = 1
                    continue
                else:
                    nreplicas += 1
                    for key, values in dict_nmolecules.items():
                        dict_nmolecules[key].append(dict_nmolecules[key][0])

                c = new_topo.add_chain()
                # Molecules
                for imol in mols_single:
                    # Just check the first atom in each molecule. It is supposed that
                    # the rest belong to the same chainID(molecule)
                    idx = list(imol)[0].index
                    name = '%s%s' % ('system', atom_to_chainid_single[idx])
                    mols_replicate.append(imol)

                # Residues
                local_ind_atom = 0
                for residue in topo_single.residues:
                    if residue.name == '':
                        name = "{0:03d}".format(ires)
                    else:
                        name = residue.name
                    dict_nresidues[name] += 1
                    list_nresidues.append(name)
                    resSeq += 1
                    r = new_topo.add_residue(name, c, resSeq, residue.segment_id)
                    # Atoms
                    for atom in residue.atoms:
                        newatom = new_topo.add_atom(atom.name, atom.element, r, serial=atom.serial)
                        newatom.serial = newatom.index + 1
                        new_xyz[0, ind_atom, :] = new_xyz[0, local_ind_atom, :] +\
                            np.dot([inx, iny, inz], trj_single.unitcell_vectors)
                        tempFactor[ind_atom] = tempFactor[local_ind_atom]
                        occupFactor[ind_atom] = occupFactor[local_ind_atom]
                        local_ind_atom += 1
                        ind_atom += 1
                    ires += 1
                # # Bonds
                for bond in topo_single.bonds:
                    a1, a2 = bond
                    a1n = new_topo.atom(a1.index+inx*nat_single+nxx*iny*nat_single+nxx*nyy*inz*nat_single)
                    a2n = new_topo.atom(a2.index+inx*nat_single+nxx*iny*nat_single+nxx*nyy*inz*nat_single)
                    new_topo.add_bond(a1n, a2n, bond.type, bond.order)

    # Check topologies ==========
    if new_topo.n_bonds != nbonds_single * nreplicas:
        mm = "Number of bonds in new topology ({}) is different to nbonds_single * nreplicas ({})".\
             format(new_topo.n_bonds, nbonds_single * nreplicas)
        print(mm) if logger_log is None else logger_log.error(mm)
        exit()

    # Save replicated pdb ==========
    if save_single:
        new_single_pdb = os.path.splitext(single_mol_pdb)[0] + "_new.pdb"
        trj_single.save_pdb(new_single_pdb)

    # Define the new box in the trajectory ==========
    new_trj = md.Trajectory(new_xyz, new_topo)
    new_trj.unitcell_lengths = np.zeros([nframes_single, 3])
    new_trj.unitcell_lengths[0, :] = [nxx * trj_single.unitcell_lengths[0, 0],
                                      nyy * trj_single.unitcell_lengths[0, 1],
                                      nzz * trj_single.unitcell_lengths[0, 2]]
    new_trj.unitcell_angles = trj_single.unitcell_angles

    # Print info of the replicated system ==========
    nows = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m1 = "\t\tMolecule replicated info ({})\n".format(nows)
    m = "\t\t" + len(m1) * "*"
    print(m1 + m) if logger_log is None else logger_log.info(m1 + m)
    nmols_replicate = new_trj.topology.find_molecules()
    m = "\t\tNumber of molecules in cell : {}\n".format(len(mols_replicate))
    m += "\t\tNumber of chains in cell    : {}\n".format(new_trj.n_chains)
    m += "\t\tNumber of residues in cell  : {}\n".format(new_trj.n_residues)
    m += "\t\tNumber of atoms in cell     : {}\n".format(new_trj.n_atoms)
    m += "\t\tSimulation box dimension    : {0:.4f} {1:.4f} {2:.4f} nm\n".\
        format(new_trj.unitcell_lengths[0, 0],
               new_trj.unitcell_lengths[0, 1],
               new_trj.unitcell_lengths[0, 2])
    m += "\t\tSimulation box angles       : {0:.2f} {1:.2f} {2:.2f} degrees\n".\
        format(new_trj.unitcell_angles[0, 0],
               new_trj.unitcell_angles[0, 1],
               new_trj.unitcell_angles[0, 2])
    m += "\t\tMass density                : {0:.4f} g/cm3\n".format(md.density(new_trj)[0]/1000)  # mdtraj gives density in kg/cm3
    m += "\t\t" + len(m1) * "*" + "\n"
    print(m) if logger_log is None else logger_log.info(m)

    # Save pdb ==========
    base, ext = os.path.splitext(single_mol_pdb)
    base = os.path.split(base)[-1]
    dirlocal = "./"
    filenamepdb_new = os.path.join(dirlocal, base + "_replicate" + ext)
    new_trj.save_pdb(filenamepdb_new, bfactors=tempFactor, )
    filenamepdb = check_and_remove_ter_labels(filenamepdb_new, occup=occupFactor, logger=logger_log)

    # Save gro =========
    # Correct number of residues > 100000
    for ires in new_trj.topology.residues:
        tmp = ires.resSeq % 100000
        ires.resSeq = tmp
    ext = ".gro"
    filenamegro_new = os.path.join(dirlocal, base + "_replicate" + ext)
    new_trj.save_gro(filenamegro_new)

    # Modify top file if exists ===========
    try:
        with open(base+".top", 'r') as f:
            lines = f.readlines()
            # Find [molecules] label index in lines
            mol_label = [s for s in lines if re.findall(r"[[ ]*molecules[ ]*]", s)]
            if len(mol_label) != 1:
                m = "The [molecules] label must appear only one time in the top file ({}).".format(base+".top")
                print(m) if logger_log is None else logger_log.error(m)
                exit()
            idx_mol_label = lines.index(mol_label[0])
            new_lines = lines[0:idx_mol_label+1]
            if nres_single/nmols_single == 1:
                # CJ2 for name_mol, n in dict_nresidues.items():
                for name_mol in list_nresidues:
                    ll = "{} {}\n".format(name_mol, 1)
                    new_lines.append(ll)
            else:
                for ireplicate in range(0, nreplicas):
                    for keys, values in dict_nmolecules.items():
                        ll = "{} {}\n".format(keys, values[ireplicate])
                        new_lines.append(ll)

        with open(base + "_replicate.top", 'w') as f:
            f.writelines(new_lines)
    except FileNotFoundError:
        pass

    if os.path.isfile("mol_noreplicate.top"):
        with open("mol_noreplicate.top", 'r') as fnr:
            lines = fnr.readlines()
            # Get name and molecules
            pos_molecules = [index for index, iline in enumerate(lines)
                             if re.match(r"^\[(\s*)molecules(\s*)\]", iline)]
            pos_molecules_new = [index for index, iline in enumerate(new_lines)
                                 if re.match(r"^\[(\s*)molecules(\s*)\]", iline)]
            with open("mol_replicate.top", "w") as ft:
                for iline in lines[:pos_molecules[0]]:
                    ft.writelines(iline)
                for iline in new_lines[pos_molecules_new[0]:]:
                    ft.writelines(iline)

    # PDB modify chain info. Mdtraj label each chain with a different name.
    # In our pdb, we want that a chain of the
    # same type have the same name.
    _modify_chain_info_pdb(filenamepdb_new, atom_to_chainid_single)

    # Map : residue_map_index[0] = [0, 1, 2, ...] Indexes of the residue 0
    residue_map_index = defaultdict(list)
    for ires in topo_single.residues:
        for iatres in ires.atoms:
            residue_map_index[ires.index].append(iatres.index)

    # Write end-to-end and backbone information
    _write_etoe_bb_info(nmols_replicate , tempFactor, occupFactor)

    # Debug
    # new_trj.save_pdb("test.pdb")
    # print(name_mol, name_residues)
    # Debug

    return trj_single, new_trj, residue_map_index, mols_replicate


# ==================================================
def _write_etoe_bb_info(nmols_replicate, beta, occup):

    fee = open("listendtoend_replicate.dat", 'w')
    fbb = open("backbone_idx_replicate.dat", 'w')
    faa = open("allatom_idx_replicate.dat", 'w')
    fee.writelines("# ich head tail\n")

    chain_to_atom = defaultdict(list)
    ich = 0
    for item in nmols_replicate:
        fbb.writelines("[mol{}]\n".format(ich))
        faa.writelines("[mol{}]\n".format(ich))
        bb_ch = []
        aa_ch = []
        for iat in item:
            idx = iat.index
            if occup[idx] == 1:
                ihead = idx
            elif occup[idx] == 2:
                itail = idx
            if beta[idx] == 0:
                bb_ch.append(idx)
            aa_ch.append(idx)

        for ibb in sorted(bb_ch):
            fbb.writelines("{}\n".format(ibb))
        for iaa in sorted(aa_ch):
            faa.writelines("{}\n".format(iaa))
        try:
            fee.writelines("{} {} {}\n".format(ich, ihead, itail))
        except Exception as e:
            continue
        ich += 1

    fee.close()
    fbb.close()


# ==================================================
def _assign_chain_to_atom(filenamepdb):

    start = ord('A')
    end = ord('Z')
    j = 1
    map_letter_to_in = defaultdict(int)
    for i in range(start, end+1):
        map_letter_to_in[chr(i)] = j
        j += 1
    map_letter_to_in[' '] = 0

    d = defaultdict()
    idx = 0
    with open(filenamepdb, 'r') as fpdb:
        lines = fpdb.readlines()
        for iline in lines:
            if iline.find("ATOM") != -1 or iline.find("HETATM") != -1:
                # idx = int(iline[6:11]) - 1
                chain_id = iline[21:22]
                d[idx] = map_letter_to_in[chain_id]
                idx += 1
    return d


# ==================================================
def _modify_chain_info_pdb(filenamepdb, atom_to_chainid_single):

    start = ord('A')
    end = ord('Z')
    j = 1
    map_in_to_letter = defaultdict(int)
    map_in_to_letter[0] = 'A'
    for i in range(start, end+1):
        map_in_to_letter[j] = chr(i)
        j += 1

    atom_to_chainid_replicate = defaultdict()
    nwlines = []
    ich = 0
    with open(filenamepdb, 'r') as fpdb:
        lines = fpdb.readlines()
        idx_local = 0
        idx = 1
        for iline in lines:
            if iline.find("ATOM") != -1 or iline.find("HETATM") != -1:
                # idx = int(iline[6:11]) - 1
                if idx > 99999:
                    magic = 100000
                    new_idx = (idx % magic) + 1
                    sentinel_a = 12
                    sentinel_b = 22
                else:
                    new_idx = idx
                    sentinel_a = 11
                    sentinel_b = 21
                try:
                    # jline = iline[0:21] + map_in_to_letter[atom_to_chainid_single[idx_local]] + iline[22:]
                    jline = iline[0:6] + "{0:5d}".format(new_idx) + iline[sentinel_a:sentinel_b] + \
                            map_in_to_letter[atom_to_chainid_single[idx_local]] + iline[sentinel_b+1:]
                    nwlines.append(jline)
                    atom_to_chainid_replicate[idx] = ich
                except KeyError:
                    idx_local = 0
                    ich += 1
                    atom_to_chainid_replicate[idx] = ich
                    # jline = iline[0:21] + map_in_to_letter[atom_to_chainid_single[idx_local]] + iline[22:]
                    jline = iline[0:6] + "{0:5d}".format(new_idx) + iline[sentinel_a:sentinel_b] + \
                            map_in_to_letter[atom_to_chainid_single[idx_local]] + iline[sentinel_b+1:]
                    nwlines.append(jline)
                idx_local += 1
                idx += 1
            elif iline.find("TER") != -1:
                # jline = iline[0:21] + map_in_to_letter[atom_to_chainid_single[idx_local - 1]] + iline[22:]
                jline = iline[0:6] + "{0:5d}".format(new_idx) + iline[sentinel_a:sentinel_b] +\
                        map_in_to_letter[atom_to_chainid_single[idx_local]] + iline[sentinel_b+1:]
                nwlines.append(jline)
                idx += 1
            else:
                nwlines.append(iline)

    with open(filenamepdb, 'w') as fpdb:
        for iline in nwlines:
            fpdb.writelines(iline)

    return atom_to_chainid_replicate


# ==================================================
def _guessing_bonds(fnameinp, topo_single, logger=None):

    # Guessing bonds from md.Universe
    try:
        m = "\n\t\t CONECT section is not present in the PDB\n"
        m += "\t\t Guessing bonds from MD.Universe"
        print(m) if logger is None else logger.info(m)
        t = mdanalysis.Universe(fnameinp, guess_bonds=True)
    # If universe cannot be created raise an error
    except Exception as e:
        m = "\n\t\t ERROR:" + str(e) + "\n"
        m += "n\t\t ERROR: Return None object from md.Universe"
        print(m) if logger is None else logger.error(m)
        exit()

    # Add bonds to the md_topology
    for ibond in t.bonds:
        idx1 = ibond.atoms[0].index
        idx2 = ibond.atoms[1].index
        iatom1 = topo_single._atoms[idx1]
        iatom2 = topo_single._atoms[idx2]
        topo_single.add_bond(iatom1, iatom2)