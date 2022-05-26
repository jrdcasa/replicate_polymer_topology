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
