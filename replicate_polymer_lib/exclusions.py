import os
import datetime


# =============================================================================
def insert_exclusions(opts, fnamepdb, residues_map_atom, logger=None):

    """

    Args:
        logger:
        residues_map_atom:
        fnamepdb:
        opts:

    Returns:

    """

    base, _ = os.path.splitext(fnamepdb)
    base = os.path.split(base)[-1]
    dirlocal = "./"

    if opts.npairs is not None:
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\tInserting exclusions in the topology file ({})\n".format(now)
        ext = ".top"
        fname = os.path.join(dirlocal, base + "_replicate" + ext)
        with open(fname, "r") as f:
            contents = f.readlines()

        idx_to_insert = contents.index("[ system ]\n")

        lines_exclusions = exclusion_gromacs(opts.npairs, residues_map_atom)

        contents.insert(idx_to_insert-1, lines_exclusions)

        with open(fname, "w") as f:
            contents = "".join(contents)
            f.write(contents)

        print(m) if logger is None else logger.info(m)


# =============================================================================
def exclusion_gromacs(npairs, residues_map_atom, logger_log=None):
    """

    Args:
        npairs:
        residues_map_atom:
        logger_log:

    Returns:

    """

    lines_exclusions = "\n[ exclusions ]\n"
    lines_exclusions += ";Number of included residues are {} neighbour(s)\n".format(npairs)
    nresidues = len(residues_map_atom)

    # For each residue
    for ires_idx, ires_values in residues_map_atom.items():

        # ires_excluded --> indexes of the residues to be excluded
        ires_idx_left = ires_idx - npairs - 1
        ires_idx_right = ires_idx + npairs + 1
        ires_excluded = []
        for i in range(0, ires_idx_left+1):
            ires_excluded.append(i)
        for i in range(ires_idx_right, nresidues):
            ires_excluded.append(i)

        # Exclusions for each atom in the excluded residues
        for iat in ires_values:
            lines_exclusions += "{} ".format(iat+1)
            icut = 0
            for jres in ires_excluded:
                for jat in residues_map_atom[jres]:
                    lines_exclusions += "{} ".format(jat+1)
                    icut += 1
                    if icut > 100:                  # Cut lines
                        lines_exclusions += "\n{} ".format(iat+1)
                        icut = 0
            lines_exclusions += "\n"

    return lines_exclusions
