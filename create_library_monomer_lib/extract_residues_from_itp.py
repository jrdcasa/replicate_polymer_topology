import re
import matplotlib.pyplot as plt
import networkx as nx
from collections import defaultdict


# =============================================================================
def extract_residues_from_itp(itpfile, wk="./", debug=False):

    """
    Extract the residues from the itp as a nx.graph

    Args:
        itpfile:
        isdraw:

    Returns:

    """

    # Initialize
    residues_dict = defaultdict(list)
    labels_graph = defaultdict(list)
    graph_itp = nx.DiGraph()

    # Open file
    is_atoms_label = False
    is_bonds_label = False
    is_moleculetype_label = False
    bonds = []
    with open(itpfile, "r") as fitp:

        lines = fitp.readlines()
        # Find the index line for r'\[\s*atoms\s*\]'
        for idx, iline in enumerate(lines):
            # Find match in the current line
            if re.findall(r'\[\s*atoms\s*]', iline):
                idx_atoms = idx
                is_atoms_label = True

            if re.findall(r'\[\s*bonds\s*]', iline):
                idx_bonds = idx
                is_bonds_label = True

            if re.findall(r'\[\s*moleculetype\s*]', iline):
                idx_moleculetype = idx
                is_moleculetype_label = True

        if not is_atoms_label:
            return None, None

        # Compile the information of the residues and then store them as attributes in the graph
        # Get info from moleculetype label
        if is_moleculetype_label:
            # Extract info
            idx = idx_moleculetype + 1
            while True:
                try:
                    iline = lines[idx]
                except IndexError:
                    break
                if re.match(r'^\[', iline) or len(iline.strip()) < 4:
                    break
                if re.match(r'^;', iline):
                    idx += 1
                    continue
                name, nrexcl = iline.split()
                idx += 1

        # Get info from atoms label
        if is_atoms_label:
            # Extract info
            idx = idx_atoms + 1
            ext = 1
            nres_current = 1
            while True:
                try:
                    iline = lines[idx]
                except IndexError:
                    break
                if re.match(r'^\[', iline) or len(iline.strip()) < 4:
                    break
                if re.match(r'^;', iline):
                    idx += 1
                    continue
                iat, nameat, nr, namer, element, cgnr, charge, mass, *_ = iline.split()
                if nres_current != int(nr):
                    nres_current = int(nr)
                    ext = 1
                label = "{0:03d}_".format(int(nr))+namer
                labels_graph[int(iat)] = label
                if element[0] == "_":
                    letter = element[1]
                else:
                    letter = element[0]
                str_at = letter+namer[0:2]+"{0:02d}".format(ext)
                node_attributes = {"labels": label, "iat": int(iat), "typeat": nameat,
                                   "nr": int(nr), "namer": namer,
                                   "element": element, "cgnr": int(cgnr),
                                   "charge": float(charge), "mass": float(mass),
                                   "nrexcl": int(nrexcl), "nameat": str_at,
                                   "label_graph": label+"_"+"{0:03d}".format(int(iat))+"_{0:s}".format(element)}
                residues_dict[int(iat)].append(node_attributes)
                graph_itp.add_node(int(iat), **node_attributes)
                idx += 1
                ext += 1

        # Get info from bonds label
        if is_bonds_label:
            # Extract info
            idx = idx_bonds + 1
            while True:
                try:
                    iline = lines[idx]
                except IndexError:
                    break
                if re.match(r'^\[', iline) or len(iline.strip()) < 4:
                    break
                if re.match(r'^;', iline):
                    idx += 1
                    continue
                try:
                    iat, jat, funct, lb0, kb = iline.split()
                except ValueError:
                    iat, jat, funct = iline.split()
                bonds.append((int(iat), int(jat)))
                idx += 1

        # Create a graph of the itp molecule
        graph_itp.add_edges_from(bonds)

        # Draw the graph
        if debug:
            pos = nx.nx_agraph.graphviz_layout(graph_itp)
            nx.draw(graph_itp, pos, with_labels=True,
                    labels=nx.get_node_attributes(graph_itp, 'label_graph'), node_size=100,
                    node_color='skyblue', font_size=4, font_color='black', font_weight='bold',
                    edge_color='gray', width=2)
            plt.show()

        return residues_dict, graph_itp
