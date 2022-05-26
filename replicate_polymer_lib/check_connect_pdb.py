import parmed as pmd
import os


# =============================================================================
def check_conect_pdb(fnamepdb):

    """
    Read the pdb file and then check if the CONECT section is present. If not the conect section is written
    and after the coordinates. The PDB file must have RESNAME fields not empty on the
    contrary a None object is returned.

    Args:
        fnamepdb (str): Path to a pdb file

    Returns:
        A string with the name of the modified pdb or None in case that parmed cannot read the pdb file

    """

    with open(fnamepdb, 'r') as f:
        lines = f.readlines()
        n = len([i for i in lines if "CONECT" in i])
        if n == 0:
            try:
                t = pmd.load_file(fnamepdb)
            except Exception as e:
                print(e)
                print("Return None object")
                return None
            graph_dict = {x: [] for x in range(0, len(t.atoms))}

            basename = os.path.splitext(os.path.split(fnamepdb)[-1])[0]
            filenamepdb = "./"+basename+"_con.pdb"

            t.save(filenamepdb, overwrite=True)
            for ibond in t.bonds:
                graph_dict[ibond.atom1.idx].append(ibond.atom2.idx)
                graph_dict[ibond.atom2.idx].append(ibond.atom1.idx)

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

            filenamepdb = fnamepdb

    return filenamepdb


# =============================================================================
def check_bonds(parmed_topo, mdtraj_topo, logger=None):

    """
    Parmed library assign bonds by comparing to standard residues templates, then based on distance criteriia for
    any atom not specified in the templates. mdtraj_topo is built using the CONECT records of a PDB file.
    The replication method can suffer of overlapped atoms and then Parmed creates new bonds. This is problematic
    when types are assigned to atoms and bonds.

    This check stops the program.

    Args:
        parmed_topo (Structure):
        mdtraj_topo (mdtraj Topology):


    """

    if len(parmed_topo.bonds) != mdtraj_topo.n_bonds:
        bl1 = []
        for ibond in mdtraj_topo.bonds:
            if ibond.atom1.index < ibond.atom2.index:
                bl1.append([ibond.atom1.index, ibond.atom2.index])
            else:
                bl1.append([ibond.atom2.index, ibond.atom1.index])
        bl2 = []
        for idx_bond in range(len(parmed_topo.bonds)):
            at1 = parmed_topo.bonds[idx_bond].atom1.idx
            at2 = parmed_topo.bonds[idx_bond].atom2.idx
            if at1 < at2:
                bl2.append([at1, at2])
            else:
                bl2.append([at2, at1])
        bond_problems = []
        for idx_bond in range(len(bl2)):
            if bl2[idx_bond] not in bl1:
                bond_problems.append(bl2[idx_bond])

        mm = "Number of bonds in topology with bond distances 'parmed' ({}) is different to real nbonds ({})\n".\
             format(len(parmed_topo.bonds), mdtraj_topo.n_bonds)
        mm += "\t\tProbably there are overlapped atoms\n"
        for ibond in bond_problems:
            mm += "\t\tBond: "+str(ibond)+" (indexes)\n"
        print(mm) if logger is None else logger.error(mm)
        exit()

