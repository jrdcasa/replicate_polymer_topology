import MDAnalysis as md
from collections import defaultdict
import os
import numpy as np
import datetime


# =============================================================================
def check_resname_pdb(fnameinp, dict_namemol=None):

    """
    Check if a pdb has name for the residues or not. It the file is not pdb return the same file
    Otherwise, use the dict_namemol to assign the namemol to each atom.
    The format for the dictionary is:

        dict_namemol = {"mol1": [0, 15], "mol2": [16, 39],...}
        The mol1 is assigned to the atom 0 to 15, the mol2 to atoms 16 to 39 and so on.

    If dict_namemol is none, the MOL name is given to all atoms.
    The new pdb is written.

    Args:
        fnameinp:
        dict_namemol:

    Returns:

    """

    ext = os.path.splitext(os.path.split(fnameinp)[-1])[-1]
    if ext != ".pdb":
        return fnameinp

    d = defaultdict()
    if dict_namemol is not None:
        for key, limits in dict_namemol.items():
            values = [i for i in range(limits[0], limits[1]+1)]
            for i in values:
                d[i] = key

    # Try to open the file
    isnewpdb = False
    with open(fnameinp, 'r') as f:
        lines = f.readlines()
        idx = 0
        newlines = []
        for iline in lines:
            if iline[0:6].find("ATOM") != -1 or iline[0:6].find("HETATM") != -1:
                if iline[17:20] == '   ':
                    isnewpdb = True
                    try:
                        newlines.append(iline[0:17]+d[idx]+iline[20:])
                    except KeyError:
                        newlines.append(iline[0:17]+"MOL"+iline[20:])
                idx += 1
            else:
                newlines.append(iline)

    if isnewpdb:
        # Write the new pdb file
        basename = os.path.splitext(os.path.split(fnameinp)[-1])[0]
        filenamepdb = "./" + basename + "_r.pdb"
        with open(filenamepdb, 'w') as f:
            for iline in newlines:
                f.writelines(iline)
        return filenamepdb
    else:
        return fnameinp

# =============================================================================
def check_conect_pdb(fnameinp, logger=None):

    """
    Read a pdb or gro file and then check if the CONECT section is present. If not the conect section is written
    after the coordinates. The PDB or GRO file must have RESNAME fields not empty on the
    contrary a None object is returned.

    Args:
        fnameinp (str): Path to a file which can be read with MDAnalysis.Universe

    Returns:
        A string with the name of the modified pdb or None in case that parmed cannot read the pdb file

    """

    bond_list_more100K = []

    # Try to open the file
    with open(fnameinp, 'r') as f:
        lines = f.readlines()
        n = len([i for i in lines if "CONECT" in i])
        natoms = len([i for i in lines if "ATOM" in i or "HETATM" in i])
        # If there is not CONECT section or the number of atoms is greater than 99999,
        # we try to build up using MDAnalysis framework
        if n == 0 or natoms > 99999:
            try:
                m = "\n\t\t CONECT section is not present in the PDB or\n"
                m += "\t\t number of atoms is >99999 (natoms={})\n".format(natoms)
                m += "\t\t Guessing bonds from MD.Universe\n"
                m += "\t\t This can be quite time-consuming.\n"
                m += "\t\t Consider using the PSF file provided by the topology " \
                     "library as an option of replicate_polymer.\n"
                now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
                m += "\t\t\t Start Guessing bonds using MDAnalysis library. ({})\n".format(now)
                print(m) if logger is None else logger.info(m)
                t = md.Universe(fnameinp, guess_bonds=True, fudge_factor=0.5)
                now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
                m += "\t\t\t End Guessing bonds using MDAnalysis library. ({})\n".format(now)
                print(m) if logger is None else logger.info(m)
            except Exception as e:
                m = "\t\t"+str(e)+"\n"
                m += "\t\tReturn None object from md.Universe"
                print(m) if logger is None else logger.info(m)
                return None, bond_list_more100K

            graph_dict = {x: [] for x in range(0, len(t.atoms))}

            basename = os.path.splitext(os.path.split(fnameinp)[-1])[0]
            if natoms < 99999:
                filenamepdb = "./"+basename+"_con.pdb"
            else:
                filenamepdb = "./" + basename + "_nocon.pdb"

            ag = t.select_atoms("all")
            ag.write(filenamepdb, frames='all')
            for ibond in t.bonds:
                graph_dict[ibond.atoms[0].ix].append(ibond.atoms[1].ix)
                graph_dict[ibond.atoms[1].ix].append(ibond.atoms[0].ix)
                bond_list_more100K.append(sorted([ibond.atoms[0].ix, ibond.atoms[1].ix]))

            if natoms < 99999:
                with open(filenamepdb, 'a') as ftmp:
                    for key, values in graph_dict.items():
                        line = "CONECT"
                        line += "{0:5d}".format(key + 1)
                        for ival in values:
                            line += "{0:5d}".format(ival + 1)
                        line += "\n"
                        ftmp.writelines(line)

            del t

        else:

            filenamepdb = fnameinp

    return filenamepdb, bond_list_more100K

# =============================================================================
def check_and_remove_ter_labels(fnameinp, occup=None, logger=None):

    """
    MDtraj library uses the TER label when several chains are presented in the pdb file.
    This function removes and renumber the atoms without these TER labels.

    Args:
        fnameinp: Name of the pdb file
        occup: List of occcup factors
        logger: a looger to report messages

    Returns:

    """

    ext = os.path.splitext(os.path.split(fnameinp)[-1])[-1]
    if ext != ".pdb":
        return fnameinp

    # Try to open the file
    isnewpdb = False
    map_oldidx_to_new_idx = defaultdict()
    with open(fnameinp, 'r') as f:
        lines = f.readlines()
        idx = 1
        newlines = []
        iline_atom = 1
        natoms = 0
        for iline in lines:
            if iline[0:6].find("ATOM") != -1 or iline[0:6].find("HETATM") != -1:
                if occup is None:
                    jline = iline[0:6]+"{0:>5d}".format(idx) + iline[11:]
                else:
                    jline = iline[0:6] + "{0:>5d}".format(idx) + iline[11:54] + "{0:6.2f}".format(occup[idx-1])+iline[60:]
                newlines.append(jline)
                map_oldidx_to_new_idx[iline_atom] = idx
                idx += 1
                iline_atom += 1
                natoms += 1
            elif iline[0:6].find("TER") != -1:
                isnewpdb = True
                iline_atom += 1
            elif iline[0:6] == "CONECT" != -1 and natoms < 100000:
                jline = "CONECT"
                pos = 6
                iiline = iline.strip() # Remove carrier return and trim spaces
                #iiline = iline[:-1]  # Remove carrier return
                while pos < len(iiline):
                    old_idx = int(iiline[pos:pos+5])
                    new_idx = map_oldidx_to_new_idx[old_idx]
                    jline += "{0:>5d}".format(new_idx)
                    pos += 5
                jline += "\n"
                newlines.append(jline)
            else:
                newlines.append(iline)

    if isnewpdb:
        # Write the new pdb file
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\n\t\tRemoving TER labels in the pdb file of the seed molecules.({})... ({})".format(fnameinp, now)
        print(m) if logger is None else logger.info(m)
        basename = os.path.splitext(os.path.split(fnameinp)[-1])[0]
        with open(fnameinp, 'w') as f:
            for iline in newlines:
                f.writelines(iline)
        return fnameinp
    else:
        return fnameinp

# # =============================================================================
# def check_bonds(parmed_topo, mdtraj_topo, logger=None):
#
#     """
#     Parmed library assign bonds by comparing to standard residues templates, then based on distance criteriia for
#     any atom not specified in the templates. mdtraj_topo is built using the CONECT records of a PDB file.
#     The replication method can suffer of overlapped atoms and then Parmed creates new bonds. This is problematic
#     when types are assigned to atoms and bonds.
#
#     This check stops the program.
#
#     Args:
#         parmed_topo (Structure):
#         mdtraj_topo (mdtraj Topology):
#
#
#     """
#
#     if len(parmed_topo.bonds) != mdtraj_topo.n_bonds:
#         bl1 = []
#         for ibond in mdtraj_topo.bonds:
#             if ibond.atom1.index < ibond.atom2.index:
#                 bl1.append([ibond.atom1.index, ibond.atom2.index])
#             else:
#                 bl1.append([ibond.atom2.index, ibond.atom1.index])
#         bl2 = []
#         for idx_bond in range(len(parmed_topo.bonds)):
#             at1 = parmed_topo.bonds[idx_bond].atom1.idx
#             at2 = parmed_topo.bonds[idx_bond].atom2.idx
#             if at1 < at2:
#                 bl2.append([at1, at2])
#             else:
#                 bl2.append([at2, at1])
#         bond_problems = []
#         for idx_bond in range(len(bl2)):
#             if bl2[idx_bond] not in bl1:
#                 bond_problems.append(bl2[idx_bond])
#
#         mm = "Number of bonds in topology with bond distances 'parmed' ({}) is different to real nbonds ({})\n".\
#              format(len(parmed_topo.bonds), mdtraj_topo.n_bonds)
#         mm += "\t\tProbably there are overlapped atoms\n"
#         for ibond in bond_problems:
#             mm += "\t\tBond: "+str(ibond)+" (indexes)\n"
#         print(mm) if logger is None else logger.error(mm)
#         exit()

