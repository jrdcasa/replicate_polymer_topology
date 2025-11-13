import networkx as nx
import os
import copy
from collections import defaultdict
import numpy as np


# =================================================================================================
def get_all_ff_for_residues(graph_itp, atomtypes, bondtype_dict,
                            angletype_dict, dihedraltype_dict, dihedralfunc_dict,
                            impropertype_dict, improper_idx_list, defaults_dict, lj14scale, iswritepairs,
                            wk="./", debug=False, log=None):

    """
    Processes and generates force field (ff) parameters for residues and links in a molecular graph.

    Parameters:
    graph_itp (networkx.Graph): The molecular graph with node attributes.
    atomtypes (dict): Dictionary of atom types.
    bondtype_dict (dict): Dictionary of bond types.
    angletype_dict (dict): Dictionary of angle types.
    dihedraltype_dict (dict): Dictionary of dihedral types.
    dihedralfunc_dict (dict): Dictionary of dihedral functions.
    impropertype_dict (dict): Dictionary of improper types.
    improper_idx_list (list): List of improper indices.
    defaults_dict (dict): Dictionary of default parameters.
    lj14scale (float): Scaling factor for 1-4 Lennard-Jones interactions.
    wk (str, optional): Working directory for output files. Defaults to "./".
    debug (bool, optional): Flag to enable debug mode. Defaults to False.
    log (logging.Logger, optional): Logger object for logging. Defaults to None.

    Returns:
    None
    """

    # Get unique labels from the graph
    unique_residues = set(label for node, label
                          in nx.get_node_attributes(graph_itp, 'labels').items())
    unique_residues = sorted(list(unique_residues))

    #  Determine combination rule based on defaults
    if defaults_dict["defaults"][1] == 2:
        combination_rule = "lorentz-berthelot"
    elif defaults_dict["defaults"][1] == 3 or defaults_dict["defaults"][1] == 1:
        combination_rule = "geometric"
    else:
        print("ERROR combination rule unknown!!!!")
        exit()
    # Determine if generating pairs
    if defaults_dict["defaults"][2] == 'no':
        is_gen_pairs = False
    elif not iswritepairs:
        is_gen_pairs = False
    else:
        is_gen_pairs = True

    # Initialize list to store subgraphs for each unique residue
    subgraphs_list = []
    for target_label in unique_residues:
        # Extract subgraph nodes for the current residue label
        subgraph_nodes = [node for node, label in nx.get_node_attributes(graph_itp, 'labels').items()
                          if label == target_label]
        subgraph = graph_itp.subgraph(subgraph_nodes)

        sg = nx.DiGraph(subgraph)

        subgraphs_list.append(nx.Graph(sg))

        # Generate force field parameters for the current residue
        generate_ff_residue(sg, atomtypes, bondtype_dict, angletype_dict,
                            dihedraltype_dict, dihedralfunc_dict, impropertype_dict, improper_idx_list,
                            defaults_dict, is_gen_pairs, combination_rule, lj14scale, iswritepairs,
                            wkdir=wk, debug=debug, log=log)

    # Generate link residues by identifying edges between different residues
    graph_edges = {tuple(sorted(t)) for t in set(graph_itp.edges())}
    tmp_subgraph_edges = set.union(*[set(subgraph.edges()) for subgraph in subgraphs_list])
    subgraph_edges = {tuple(sorted(t)) for t in tmp_subgraph_edges}
    edges_links = sorted(graph_edges - subgraph_edges)

    # Generate force field parameters for the link residues
    generate_links_residue(graph_itp, edges_links, atomtypes, bondtype_dict,
                           angletype_dict, dihedraltype_dict, dihedralfunc_dict,
                           impropertype_dict, improper_idx_list, defaults_dict,
                           is_gen_pairs, combination_rule, lj14scale, iswritepairs,
                           wk=wk, log=log)


# =================================================================================================
def generate_ff_residue(sg, atomtypes, bondtype_dict, angletype_dict,
                        dihedraltype_dict, dihedraltype_function_dict, impropertype_dict,
                        improper_idx_list, defaults_dict,
                        is_gen_pairs, combination_rule, lj14scale, iswritepairs,
                        wkdir="./", debug=False, log=None):

    """
    Generate the terms of the ff for polyply

    Args:
        sg:
        atomtypes
        bondtype_dict:
        angletype_dict:
        dihedraltype_dict:
        dihedraltype_function_dict:
        impropertype_dict:
        improper_idx_list:
        defaults_dict:
        is_gen_pairs:
        combination_rule:
        lj14scale
        wkdir:
        debug:
        log:

    Returns:

    """

    # Get the root node
    root_node = next(nx.topological_sort(sg))

    labels = nx.get_node_attributes(sg, 'namer')
    nrexcl = nx.get_node_attributes(sg, 'nrexcl')[root_node]
    pattern = labels[root_node]

    polypy_fffile = os.path.join(wkdir, pattern+".ff")
    if not os.path.isfile(polypy_fffile):

        with open(polypy_fffile, "w") as fp:

            lines = "[ moleculetype ]\n"
            lines += "; Name      nrexcl\n"
            lines += "{0:s}         {1:d}\n".format(pattern, nrexcl)

            # Write atoms
            lines += write_polyply_atoms(sg, root_node)

            # Get bonds
            lines += write_polyply_bonds(sg, root_node, bondtype_dict, wkdir, log)

            # Get angles
            lines += write_polyply_angles(sg, root_node, angletype_dict, wkdir, log)

            # Get dihedrals
            lines += write_polyply_dihedrals(sg, root_node, dihedraltype_dict, dihedraltype_function_dict, log)

            # Get impropers
            lines += write_polyply_impropers(sg, root_node, impropertype_dict, improper_idx_list, wkdir, log)

            # Get pairs
            if iswritepairs:
                lines += write_polyply_pairs(sg, root_node, atomtypes, defaults_dict,
                                             is_gen_pairs, combination_rule, lj14scale)

            # Write the lines
            fp.writelines(lines)

    # Debug =================================================================================
    if debug:
        import matplotlib.pyplot as plt
        pos = nx.nx_agraph.graphviz_layout(sg)
        nx.draw(sg, pos, with_labels=True, labels=nx.get_node_attributes(sg, 'label_graph'), node_size=700,
                node_color='skyblue', font_size=6, font_color='black',
                font_weight='bold', edge_color='gray', width=2)
        plt.show()
    # Debug =================================================================================


# =======================================================================================
def find_all_simple_paths(graph, length):

    def find_paths(graphf, start, target_length):

        stack = [(start, [start])]
        pathsf = []

        while stack:
            current_node, current_path = stack.pop()

            if len(current_path) == target_length:
                pathsf.append(current_path)
                continue

            for neighbor in graphf.neighbors(current_node):
                if neighbor not in current_path:
                    stack.append((neighbor, current_path + [neighbor]))

        return pathsf

    # Convert the directed graph to an undirected graph
    undirected_graph = nx.Graph(graph)
    all_paths_of_length = []
    for node in graph.nodes:
        paths = find_paths(undirected_graph, node, length)
        all_paths_of_length.extend(paths)

    # Remove repeats paths
    for ipath in all_paths_of_length:
        ipath_r = ipath[::-1]
        if ipath_r in all_paths_of_length:
            all_paths_of_length.remove(ipath_r)

    return all_paths_of_length


# =======================================================================================
def find_node_simple_paths(graph, node_find, length):

    def find_paths(graphf, start, target_length):

        stack = [(start, [start])]
        pathsf = []

        while stack:
            current_node, current_path = stack.pop()

            if len(current_path) == target_length:
                pathsf.append(current_path)
                continue

            for neighbor in graphf.neighbors(current_node):
                if neighbor not in current_path:
                    stack.append((neighbor, current_path + [neighbor]))

        return pathsf

    # Convert the directed graph to an undirected graph
    undirected_graph = nx.Graph(graph)
    all_node_paths_of_length = []
    paths = find_paths(undirected_graph, node_find, length)
    all_node_paths_of_length.extend(paths)

    # Remove repeats paths
    for ipath in all_node_paths_of_length:
        ipath_r = ipath[::-1]
        if ipath_r in all_node_paths_of_length:
            all_node_paths_of_length.remove(ipath_r)

    return all_node_paths_of_length


# =======================================================================================
def write_polyply_atoms(sg, root_node):

    # This not work for system with branches. Gives the order ot atoms incorrect
    bfs_tree = list(nx.bfs_tree(sg, root_node))
    # Order the number of the nodes. The nodes have been already order following the Topology rules.
    bfs_tree = sorted(bfs_tree)

    idx = 1
    atom_lines = "[ atoms ]\n"
    atom_lines += ";   nr       type  resnr residue  atom   cgnr     charge       mass\n"
    sum_charge = 0.0
    # Loop through nodes in BFS order and get attributes
    visited_nodes = set()
    for node in bfs_tree:
        if node not in visited_nodes:
            visited_nodes.add(node)
            node_attributes = sg.nodes[node]
            sum_charge += float(node_attributes['charge'])
            atom_lines += "{0:>6d}{1:>12s}{2:8d}{3:>8s}{4:>8s}{5:>6d}{6:12.4f}{7:12.4f}   ; {8:7.4f}\n" \
                .format(idx, node_attributes['typeat'], 1, node_attributes['namer'],
                        node_attributes['nameat'], idx,
                        node_attributes['charge'], node_attributes['mass'], sum_charge)
            idx += 1
            # Debug =================================================================================
            # print(f"Node {node} attributes:", node_attributes)
            # Debug =================================================================================

    return atom_lines


# =======================================================================================
def write_polyply_bonds(sg, root_node, bondtype_dict, wkdir, log):

    bfs_edges = list(nx.bfs_edges(sg, root_node))

    bond_lines = "[ bonds ]\n"
    bond_lines += "; atom1 atom2   b0    kb\n"
    for ibond in bfs_edges:
        i, j = ibond
        node_attributes_i = sg.nodes[i]
        node_attributes_j = sg.nodes[j]
        itype = node_attributes_i['typeat']
        jtype = node_attributes_j['typeat']
        labelf = itype + "-" + jtype
        labelr = jtype + "-" + itype
        if labelf in bondtype_dict.keys():
            label = labelf
        elif labelr in bondtype_dict.keys():
            label = labelr
        else:
            label = ""
            m = "\n\t\t ERROR: Bond type {}.\n".format(labelf)
            m += "\t\t ERROR: or Bond type {}.\n".format(labelr)
            m += "\t\t ERROR: don't have kbond parameters\n"
            m += "\t\t ERROR: Working directory: {}\n".format(wkdir)
            print(m) if log is None else log.error(m)
            exit()
        inew = i - root_node + 1
        jnew = j - root_node + 1
        function_bond = bondtype_dict[label][0][-1]
        if function_bond != 3:
            bond_lines += "{0:>6d}{1:>6d}{2:>6d}{3:>12.4f}{4:>12.1f}\n".format(inew, jnew, function_bond,
                                                                               bondtype_dict[label][0][0],
                                                                               bondtype_dict[label][0][1])
        else:
            m = "\n\t\tERROR: Morse potential is not implemented.\n"
            m += "\t\t\tget_all_ff_for_residues.py. write_polyply_bonds()\n"
            print(m) if log is None else log.error(m)
            exit()
    return bond_lines


# =======================================================================================
def write_polyply_angles(sg, root_node, angletype_dict, wkdir, log):

    angles = find_all_simple_paths(sg, 3)
    angle_lines = "[ angles ]\n"
    angle_lines += "; ai    aj    ak funct  theta kbend\n"
    for iangle in angles:
        i, j, k = iangle
        node_attributes_i = sg.nodes[i]
        node_attributes_j = sg.nodes[j]
        node_attributes_k = sg.nodes[k]
        itype = node_attributes_i['typeat']
        jtype = node_attributes_j['typeat']
        ktype = node_attributes_k['typeat']
        labelf = itype + "-" + jtype + "-" + ktype
        labelr = ktype + "-" + jtype + "-" + itype
        if labelf in angletype_dict.keys():
            label = labelf
        elif labelr in angletype_dict.keys():
            label = labelr
        else:
            label = ""
            m = "\n\t\t ERROR: Angle type {}.\n".format(labelf)
            m += "\t\t ERROR: or Angle type {}.\n".format(labelr)
            m += "\t\t ERROR: don't have kbend parameters\n"
            m += "\t\t ERROR: Working directory: {}\n".format(wkdir)
            print(m) if log is None else log.error(m)
            exit()
        inew = i - root_node + 1
        jnew = j - root_node + 1
        knew = k - root_node + 1
        # For Harmonic and Cosine angle forms (type 1 and 2)
        function_angle = angletype_dict[label][0][-1]
        if function_angle == 1 or function_angle ==2:
            angle_lines += "{0:>6d}{1:>6d}{2:>6d}{3:>6d}{4:>12.4f}{5:>12.1f}\n".format(inew, jnew, knew,
                                                                                       function_angle,
                                                                                       angletype_dict[label][0][0],
                                                                                       angletype_dict[label][0][1])
        elif function_angle == 5:  # Urey-Bradley form (type 5)
            angle_lines += "{0:>6d}{1:>6d}{2:>6d}{3:>6d} {4:>12.4f} {5:>12.1f} {6:>12.4f} {7:>12.1f}\n".\
                format(inew, jnew, knew,
                       function_angle ,
                       angletype_dict[label][0][0],
                       angletype_dict[label][0][1],
                       angletype_dict[label][0][2],
                       angletype_dict[label][0][3])
        else:
            m = "\n\t\tERROR: Angle potential {} is not implemented.\n".format(function_angle)
            m += "\t\t\tget_all_ff_for_residues.py. write_polyply_angles()\n"
            print(m) if log is None else log.error(m)
            exit()

    return angle_lines


# =======================================================================================
def write_polyply_dihedrals(sg, root_node, dihedraltype_dict, dihedraltype_function_dict, log):

    dihedrals = find_all_simple_paths(sg, 4)
    dihedral_lines = "[ dihedrals ]\n"
    dihedral_lines += "; propers\n"
    dihedrals_notfound = defaultdict()
    ndihedrals_nottyped = 0
    ndihedrals_nottyped_noth = 0
    for idih in dihedrals:
        i, j, k, m = idih
        node_attributes_i = sg.nodes[i]
        node_attributes_j = sg.nodes[j]
        node_attributes_k = sg.nodes[k]
        node_attributes_m = sg.nodes[m]
        itype = node_attributes_i['typeat']
        jtype = node_attributes_j['typeat']
        ktype = node_attributes_k['typeat']
        mtype = node_attributes_m['typeat']
        labelf = itype + "-" + jtype + "-" + ktype + "-" + mtype
        labelr = mtype + "-" + ktype + "-" + jtype + "-" + itype
        if labelf in dihedraltype_dict.keys():
            label = labelf
        elif labelr in dihedraltype_dict.keys():
            label = labelr
        else:
            ndihedrals_nottyped += 1
            if node_attributes_i['element'] != 'H' and \
                    node_attributes_m['element'] != 'H':
                ndihedrals_nottyped_noth += 1
            if labelf not in dihedrals_notfound.keys() and \
                    labelr not in dihedrals_notfound.keys():
                dihedrals_notfound[labelf] = labelr
            continue
        inew = i - root_node + 1
        jnew = j - root_node + 1
        knew = k - root_node + 1
        mnew = m - root_node + 1
        dihedral_label_lines = "{0:>6d}{1:>6d}{2:>6d}{3:>6d}{4:>6d}".format(inew, jnew, knew, mnew,
                                                                            dihedraltype_function_dict[label])
        if dihedraltype_function_dict[label] == 3:
            # Ryckaert-Bellemans dihedral
            ll = dihedral_label_lines
            for i in dihedraltype_dict[label][0][:-1]:
                ll += "{0:>12.3f}".format(i)
            dihedral_lines += ll+"\n"
        elif dihedraltype_function_dict[label] == 9:
            for idx, lparm in enumerate(dihedraltype_dict[label]):
                ll = dihedral_label_lines + " {0:>12.1f}   {1:>12.3f}   {2:>2d}\n".format(lparm[0], lparm[1], lparm[2])
                dihedral_lines += ll
        else:
            m = "\n\t\t ERROR: Dihedral function {} is not " \
                "yet implemented in get_all_for_residues.py.".format(dihedraltype_function_dict[label])
            print(m) if log is None else log.error(m)

    m = "\t\t ======== Dihedrals in fragments ========"
    print(m) if log is None else log.info(m)
    for key, values in dihedrals_notfound.items():
        m = "\t\t WARNING: Dihedral type {} (or {}) cannot be found.".format(key, values)
        print(m) if log is None else log.info(m)
    m = "\t\t Number of untyped dihedrals                                             : {}\n" \
        .format(ndihedrals_nottyped)
    m += "\t\t Number of untyped dihedrals that involve non-heavy atoms at the termini : {}\n" \
        .format(ndihedrals_nottyped_noth)
    print(m) if log is None else log.info(m)

    return dihedral_lines


# =======================================================================================
def write_polyply_impropers(sg, root_node, impropertype_dict, improper_idx_list, wkdir, log):

    """
    Write the impropers found in the itp used as input.
    Args:
        sg:
        root_node:
        impropertype_dict:
        improper_idx_list:
        wkdir:
        log:

    Returns:

    """

    improper_lines = "; impropers\n"
    for item_improper in improper_idx_list:
        found_improper = True
        label = ""
        for idx, iatom in enumerate(item_improper):
            iatom_local = iatom - (root_node - 1)
            if iatom not in sg.nodes():
                found_improper = False
                break
            else:
                label += sg.nodes[iatom]['typeat']+"-"
        if found_improper:
            label = label[0:-1]
            improper_lines += "    ".join(map(str, item_improper))
            param = impropertype_dict[label][0]
            type_improper = param[-1]
            if type_improper == 4:
                improper_lines += "{0:4d} {1:7.3f} {2:7.3f} {3:4d}".format(type_improper, param[0], param[1],
                                                                           int(param[2]))
            elif type_improper == 2:
                improper_lines += "{0:4d} {1:7.3f} {2:7.3f}".format(type_improper, param[0], param[1])

            improper_lines += "\n"

    return improper_lines


# =======================================================================================
def write_polyply_pairs(sg, root_node, atomtypes, defaults_dict, is_gen_pairs, combination_rule, lj14scale):

    # Write pairs
    atoms = list(sg.nodes())
    pair_lines = "[ pairs ]\n"
    for iat in atoms:
        iatnew = iat - root_node + 1
        # pairs = find_all_simple_paths(sg, 4)
        pairs = find_node_simple_paths(sg, iat, 4)
        for ipair in pairs:
            jat = ipair[-1]
            jatnew = jat - root_node + 1
            if iatnew > jatnew:
                continue
            if is_gen_pairs:
                pair_lines += "{0:6d}{1:6d}{2:6d}\n".format(iatnew, jatnew, 1)
            else:
                iatype = nx.get_node_attributes(sg, 'typeat')[iat]
                jatype = nx.get_node_attributes(sg, 'typeat')[jat]
                epsilon_ij = lj14scale*np.sqrt(atomtypes[iatype][-1]*atomtypes[jatype][-1])
                if combination_rule == "lorentz-berthelot":
                    sigma_ij = np.sqrt(atomtypes[iatype][-2] * atomtypes[jatype][-2])
                elif combination_rule == "geometric":
                    sigma_ij = 0.5 * (atomtypes[iatype][-2] + atomtypes[jatype][-2])
                pair_lines += "{0:6d}{1:6d}{2:6d} {3:12.9f} {4:12.9f}\n".format(iatnew, jatnew,
                                                                                defaults_dict["defaults"][0],
                                                                                sigma_ij, epsilon_ij)

    return pair_lines


# =======================================================================================
def generate_links_residue(graph_itp, edges_links, atomtypes, bondtype_dict, angletype_dict,
                           dihedraltype_dict, dihedralfunc_dict,
                           impropertype_dict, improper_idx_list,
                           defaults_dict,
                           is_gen_pairs, combination_rule, lj14scale, iswritepairs,
                           wk="./", log=None):

    graph_itp_undirected = nx.Graph(graph_itp)

    for iedge in edges_links:
        idx_tail, idx_head = iedge
        res_head = nx.get_node_attributes(graph_itp_undirected, 'labels')[idx_head].split("_")[-1]
        res_tail = nx.get_node_attributes(graph_itp_undirected, 'labels')[idx_tail].split("_")[-1]
        resname_tail = nx.get_node_attributes(graph_itp_undirected, 'labels')[idx_tail]
        resname_head = nx.get_node_attributes(graph_itp_undirected, 'labels')[idx_head]
        type_head = nx.get_node_attributes(graph_itp_undirected, 'typeat')[idx_head]
        type_tail = nx.get_node_attributes(graph_itp_undirected, 'typeat')[idx_tail]
        name_head = nx.get_node_attributes(graph_itp_undirected, 'nameat')[idx_head]
        name_tail = nx.get_node_attributes(graph_itp_undirected, 'nameat')[idx_tail]
        filename = os.path.join(wk, res_tail+"_"+res_head+"_links.ff")

        if not os.path.isfile(filename):

            with open(filename, 'w') as flink:
                m = "\t\t Creating the link file: {} ".format(filename)
                print(m) if log is None else log.info(m)
                lines = ";-------------------------------------\n"
                lines += ";             {}                      \n".format(res_tail+"-"+res_head+" LINK")
                lines += ";-------------------------------------\n"
                lines += "[link]\n"
                lines_term = ";-------------------------------------\n"
                lines_term += ";             {}                      \n".format(res_tail+"-"+res_head+" LINK-TERMS")
                lines_term += ";-------------------------------------\n"
                lines_term += "[link]\n"

                if res_head.upper() == res_tail.upper():
                    lines += 'resname "{}"\n'.format(res_head)
                    lines_term += 'resname "{}"\n'.format(res_head)
                else:
                    lines += 'resname "{}|{}"\n'.format(res_tail, res_head)
                    lines_term += 'resname "{}|{}"\n'.format(res_tail, res_head)

                # Bonds
                lines += '[ bonds ]\n'
                labelf = type_head + "-" + type_tail
                labelr = type_tail + "-" + type_head
                if labelf in bondtype_dict.keys():
                    label = labelf
                elif labelr in bondtype_dict.keys():
                    label = labelr
                else:
                    mm = "\n\t\t ERROR: Bond type {}.\n".format(labelf)
                    mm += "\t\t ERROR: or Bond type {}.\n".format(labelr)
                    mm += "\t\t ERROR: don't have kbond parameters\n"
                    mm += "\t\t ERROR: Working directory: {}\n".format(wk)
                    print(mm) if log is None else log.error(mm)
                    exit()
                lines += '+{0:<6s}{1:6s}{2:>6d}{3:>12.4f}{4:>12.1f}  {{"group": "connection"}}\n'.\
                    format(name_head, name_tail, bondtype_dict[label][0][-1],
                           bondtype_dict[label][0][0], bondtype_dict[label][0][1])

                # Angles ======================
                lines += '[ angles ]\n'
                # Neighbours of Tail atom
                for ilink_atom in [idx_head, idx_tail]:
                    for ipair in find_node_simple_paths(graph_itp_undirected, ilink_atom, 2):
                        if sorted(ipair) == sorted(iedge):
                            continue
                        ktype = nx.get_node_attributes(graph_itp_undirected, 'typeat')[ipair[1]]
                        kname = nx.get_node_attributes(graph_itp_undirected, 'nameat')[ipair[1]]
                        if ilink_atom == idx_tail:
                            labelf = type_head + "-" + type_tail + "-" + ktype
                            labelr = ktype + "-" + type_tail + "-" + type_head
                        else:
                            labelf = type_tail + "-" + type_head + "-" + ktype
                            labelr = ktype + "-" + type_head + "-" + type_tail
                        if labelf in angletype_dict.keys():
                            label = labelf
                        elif labelr in angletype_dict.keys():
                            label = labelr
                        else:
                            mm = "\n\t\t ERROR: Angle type {}.\n".format(labelf)
                            mm += "\t\t ERROR: or Angle type {}.\n".format(labelr)
                            mm += "\t\t ERROR: don't have kbond parameters\n"
                            mm += "\t\t ERROR: Working directory: {}\n".format(wk)
                            print(mm) if log is None else log.error(mm)
                            exit()
                        func_angle = angletype_dict[label][0][-1]
                        if func_angle == 1 or func_angle == 2:
                            if ilink_atom == idx_tail:
                                lines += '{0:<8s} {1:<8s}+{2:<8s} {3:<6d} {4:>12.4f} {5:>12.1f}   ' \
                                         '   {{"group": "connection"}}\n'.\
                                    format(kname, name_tail, name_head, func_angle,
                                           angletype_dict[label][0][0], angletype_dict[label][0][1])
                            else:
                                lines += '{0:<8s}+{1:<8s}+{2:<8s} {3:<6d} {4:>12.4f} {5:>12.1f}   ' \
                                         '   {{"group": "connection"}}\n'.\
                                    format(name_tail, name_head, kname, func_angle,
                                           angletype_dict[label][0][0], angletype_dict[label][0][1])
                        elif func_angle == 5:
                            if ilink_atom == idx_tail:
                                lines += '{0:<8s} {1:<8s}+{2:<8s} {3:<6d} {4:>12.4f} {5:>12.1f} {6:>12.4f} {7:>12.1f}' \
                                         '   {{"group": "connection"}}\n'.\
                                    format(kname, name_tail, name_head, func_angle,
                                           angletype_dict[label][0][0], angletype_dict[label][0][1],
                                           angletype_dict[label][0][2], angletype_dict[label][0][3])
                            else:
                                lines += '{0:<8s}+{1:<8s}+{2:<8s} {3:<6d} {4:>12.4f} {5:>12.1f} {6:>12.4f} {7:>12.1f} '\
                                         '   {{"group": "connection"}}\n'.\
                                    format(name_tail, name_head, kname, func_angle,
                                           angletype_dict[label][0][0], angletype_dict[label][0][1],
                                           angletype_dict[label][0][2], angletype_dict[label][0][3])

                # Dihedrals ======================
                lines += '[ dihedrals ]\n ; propers\n'
                lines_term += '[ dihedrals ]\n'
                link_dihedrals = list()
                mask_link_dihedral = list()
                link_pairs = list()
                mask_link_pairs = list()
                dihedrals_notfound = defaultdict()
                ndihedrals_nottyped = 0
                ndihedrals_nottyped_noth = 0
                idx_tail, idx_head = iedge
                subdihs = find_node_simple_paths(graph_itp_undirected, idx_tail, 4)
                for item in subdihs:
                    resname_nexttail = nx.get_node_attributes(graph_itp_undirected, 'labels')[item[1]]
                    if resname_tail == resname_nexttail:
                        continue
                    link_dihedrals.append(item)
                    tail_order = int(nx.get_node_attributes(graph_itp_undirected, 'labels')[idx_tail].split("_")[0])
                    tmp_list = [" "]
                    for e in item[1:]:
                        next_order = int(nx.get_node_attributes(graph_itp_undirected, 'labels')[e].split("_")[0])
                        diff = next_order - tail_order
                        if next_order > tail_order:
                            tmp_list.append(diff*"+")
                        elif next_order < tail_order:
                            tmp_list.append(diff*"-")
                        else:
                            tmp_list.append(" ")
                    mask_link_dihedral.append(tmp_list)
                    link_pairs.append([item[0], item[3]])
                    mask_link_pairs.append([" ", diff*"+"])

                subdihs = find_node_simple_paths(graph_itp_undirected, idx_tail, 3)
                for item in subdihs:
                    resname_nexttail = nx.get_node_attributes(graph_itp_undirected, 'labels')[item[1]]
                    if resname_tail == resname_nexttail:
                        continue
                    graph_itp_undirected.neighbors(idx_tail)
                    for ineigh in graph_itp_undirected.neighbors(idx_tail):
                        if ineigh == idx_head:
                            continue
                        item.insert(0, ineigh)
                        item_copy = copy.deepcopy(item)
                        link_dihedrals.append(item_copy)
                        tail_order = int(nx.get_node_attributes(graph_itp_undirected, 'labels')[idx_tail].split("_")[0])
                        tmp_list = [" "]
                        for e in item[1:]:
                            next_order = int(nx.get_node_attributes(graph_itp_undirected, 'labels')[e].split("_")[0])
                            diff = next_order - tail_order
                            if next_order > tail_order:
                                tmp_list.append(diff * "+")
                            elif next_order < tail_order:
                                tmp_list.append(diff * "-")
                            else:
                                tmp_list.append(" ")
                        mask_link_dihedral.append(tmp_list)
                        link_pairs.append([item[0], item[3]])
                        mask_link_pairs.append([" ", diff * "+"])
                        item.pop(0)

                nodes = find_node_simple_paths(graph_itp_undirected, idx_tail, 3)
                subdihs = list()
                for item in nodes:
                    ires_item0 = int(nx.get_node_attributes(graph_itp_undirected, 'labels')[item[0]].split("_")[0])
                    ires_item2 = int(nx.get_node_attributes(graph_itp_undirected, 'labels')[item[-1]].split("_")[0])
                    if ires_item0 != ires_item2:
                        continue
                    tmp_subdihs = find_node_simple_paths(graph_itp_undirected, item[-1], 4)
                    for jtem in tmp_subdihs:
                        jres_item0 = int(nx.get_node_attributes(graph_itp_undirected, 'labels')[jtem[0]].split("_")[0])
                        jres_item2 = int(nx.get_node_attributes(graph_itp_undirected, 'labels')[jtem[-1]].split("_")[0])
                        if jres_item2 <= jres_item0:
                            continue
                        link_dihedrals.append(jtem)
                        mask_link_dihedral.append([" ", " ", " ", "+"])
                        link_pairs.append([jtem[0], jtem[3]])
                        mask_link_pairs.append([" ", "+"])

                # Write all dihedrals
                for jdx, item in enumerate(link_dihedrals):
                    itype = nx.get_node_attributes(graph_itp_undirected, 'typeat')[item[0]]
                    iname = nx.get_node_attributes(graph_itp_undirected, 'nameat')[item[0]]
                    jtype = nx.get_node_attributes(graph_itp_undirected, 'typeat')[item[1]]
                    jname = nx.get_node_attributes(graph_itp_undirected, 'nameat')[item[1]]
                    ktype = nx.get_node_attributes(graph_itp_undirected, 'typeat')[item[2]]
                    kname = nx.get_node_attributes(graph_itp_undirected, 'nameat')[item[2]]
                    mtype = nx.get_node_attributes(graph_itp_undirected, 'typeat')[item[3]]
                    mname = nx.get_node_attributes(graph_itp_undirected, 'nameat')[item[3]]

                    labelf = itype + "-" + jtype + "-" + ktype + "-" + mtype
                    labelr = mtype + "-" + ktype + "-" + jtype + "-" + itype

                    if labelf in dihedraltype_dict.keys():
                        label = labelf
                    elif labelr in dihedraltype_dict.keys():
                        label = labelr
                    else:
                        ndihedrals_nottyped += 1
                        if nx.get_node_attributes(graph_itp_undirected, 'element')[item[0]] != 'H' and \
                                nx.get_node_attributes(graph_itp_undirected, 'element')[item[3]] != 'H':
                            ndihedrals_nottyped_noth += 1
                        if labelf not in dihedrals_notfound.keys() and \
                                labelr not in dihedrals_notfound.keys():
                            dihedrals_notfound[labelf] = labelr
                        continue

                    lines_tmp = "{0:1s}{1:<6s}{2:1s}{3:<6s}{4:1s}{5:<6s}{6:1s}{7:<6s}{8:>6d} ". \
                        format(mask_link_dihedral[jdx][0], iname,
                               mask_link_dihedral[jdx][1], jname,
                               mask_link_dihedral[jdx][2], kname,
                               mask_link_dihedral[jdx][3], mname,
                               dihedralfunc_dict[label])
                    header = lines_tmp

                    if dihedralfunc_dict[label] == 3:
                        for i in dihedraltype_dict[label][0][:-1]:
                            lines_tmp += "{0:>12.3f}".format(i)
                        lines_tmp += ' {"group": "connection"}\n'
                    elif dihedralfunc_dict[label] == 9:
                        for idx, i in enumerate(dihedraltype_dict[label]):
                            if idx == 0:
                                lines_tmp += " {0:>12.1f}   {1:>12.3f}   {2:>2d}\n".format(i[0], i[1], i[2])
                            else:
                                lines_tmp += header + " {0:>12.1f}   {1:>12.3f}   {2:>2d}\n".format(i[0], i[1], i[2])
                    if lines_tmp.count("++") > 0:
                        lines_term += lines_tmp
                    else:
                        lines += lines_tmp

                # # Impropers ======================
                lines += '; impropers\n'
                link_impropers = list()
                lines_tmp = ""
                for improper_idx in improper_idx_list:
                    improper_nameat = [nx.get_node_attributes(graph_itp_undirected, 'nameat')[i] for i in improper_idx]
                    improper_nr = [nx.get_node_attributes(graph_itp_undirected, 'nr')[i] for i in improper_idx]
                    improper_resnam = [nx.get_node_attributes(graph_itp_undirected, 'labels')[i] for i in improper_idx]
                    improper_typeat = [nx.get_node_attributes(graph_itp_undirected, 'typeat')[i] for i in improper_idx]
                    nodes = [node for node, label in nx.get_node_attributes(graph_itp_undirected, 'labels').items()
                             if label == resname_tail or label == resname_head]

                    if len(set(improper_nr)) > 1 and len(set(improper_resnam)) > 1 \
                            and all(i in nodes for i in improper_idx):
                        max_ir = max(improper_nr)
                        min_ir = min(improper_nr)
                        delta = max_ir - min_ir
                        for idx, ir in enumerate(improper_nr):
                            if ir == max_ir:
                                for j in range(0, delta):
                                    lines_tmp += "+"
                            lines_tmp += "{} ".format(improper_nameat[idx])
                        label_improper = "{0:s}-{1:s}-{2:s}-{3:s}".\
                            format(improper_typeat[0], improper_typeat[1], improper_typeat[2], improper_typeat[3])

                        improper_func = impropertype_dict[label_improper][0][-1]
                        if improper_func == 2:
                            lines_tmp += "{0:4d} {1:7.3f} {2:7.3f} ".\
                                format(improper_func,
                                       impropertype_dict[label_improper][0][0],
                                       impropertype_dict[label_improper][0][1])
                        else:
                            lines_tmp += "{0:4d} {1:7.3f} {2:7.3f} {3:4d} ".\
                                format(improper_func,
                                       impropertype_dict[label_improper][0][0],
                                       impropertype_dict[label_improper][0][1],
                                       int(impropertype_dict[label_improper][0][2]))
                        lines_tmp += '{"group": "connection"} \n'
                if lines_tmp.count("++") > 0:
                    lines_term += lines_tmp
                else:
                    lines += lines_tmp

                # Write pairs
                if iswritepairs:
                    lines += '[ pairs ]\n'
                    lines_term += '[ pairs ]\n'
                    for jdx, item in enumerate(link_pairs):
                        iname = nx.get_node_attributes(graph_itp_undirected, 'nameat')[item[0]]
                        jname = nx.get_node_attributes(graph_itp_undirected, 'nameat')[item[1]]
                        if is_gen_pairs:
                            lines_tmp = '{0:1s}{1:<6s}{2:1s}{3:<6s} 1 {{"group": "connection"}}\n'. \
                                format(mask_link_pairs[jdx][0], iname,
                                       mask_link_pairs[jdx][1], jname)
                        else:
                            iatype = nx.get_node_attributes(graph_itp, 'typeat')[item[0]]
                            jatype = nx.get_node_attributes(graph_itp, 'typeat')[item[1]]
                            epsilon_ij = lj14scale*np.sqrt(atomtypes[iatype][-1]*atomtypes[jatype][-1])
                            if combination_rule == "lorentz-berthelot":
                                sigma_ij = np.sqrt(atomtypes[iatype][-2] * atomtypes[jatype][-2])
                            elif combination_rule == "geometric":
                                sigma_ij = 0.5 * (atomtypes[iatype][-2] + atomtypes[jatype][-2])

                            lines_tmp = '{0:1s}{1:<6s}{2:1s}{3:<6s} {4:6d} {5:12.9f} {6:12.9f}' \
                                        ' {{"group": "connection"}}\n'. \
                                format(mask_link_pairs[jdx][0], iname,
                                       mask_link_pairs[jdx][1], jname, defaults_dict["defaults"][0], sigma_ij, epsilon_ij )

                        # lines_tmp = '{0:1s}{1:<6s}{2:1s}{3:<6s} 1 {{"group": "connection"}}\n'. \
                        #     format(mask_link_pairs[jdx][0], iname,
                        #            mask_link_pairs[jdx][1], jname)
                        if lines_tmp.count("++") > 0:
                            lines_term += lines_tmp
                        else:
                            lines += lines_tmp

                mm = "\t\t ======== Dihedrals in link ========"
                print(mm) if log is None else log.info(mm)
                for key, values in dihedrals_notfound.items():
                    mm = "\t\t WARNING: Dihedral type {} (or {}) cannot be found.".format(key, values)
                    print(mm) if log is None else log.info(mm)
                mm = "\t\t Number of untyped dihedrals                                             : {}\n" \
                    .format(ndihedrals_nottyped)
                mm += "\t\t Number of untyped dihedrals that involve non-heavy atoms at the termini : {}\n" \
                    .format(ndihedrals_nottyped_noth)
                print(mm) if log is None else log.info(mm)

                flink.writelines(lines)
                flink.writelines(lines_term)

    m = "\n"
    print(m) if log is None else log.info(m)
