from topology import Topology
from topology.atomic_data import atomic_number
from replicate_polymer_lib.forcefield.Smarts_Graph import SMARTSGraph
from replicate_polymer_lib.forcefield.topology_graph import TopologyGraph
import ele
from ele.exceptions import ElementError


def _convert_topology_to_topologygraph(topology):

    top_graph = TopologyGraph()

    # Atoms =========
    for iatom_idx in topology._graphdict.keys():
        # For no-real atoms, rather interaction sites (i.e: CG or UA models)
        if topology._names[iatom_idx].startswith("_"):
            top_graph.add_atom(name=topology._names[iatom_idx],
                               index=iatom_idx,
                               atomic_number=None,
                               element=None)
        else:
            top_graph.add_atom(name=topology._names[iatom_idx],
                               index=iatom_idx,
                               atomic_number=atomic_number[topology.elements[iatom_idx]],
                               element=topology.elements[iatom_idx])

    # Bonds =========
    for ibond in topology._bonds:
        ll = list(ibond)
        idx = ll[0]
        jdx = ll[1]
        if jdx > idx:
            top_graph.add_bond(jdx, idx)
        else:
            top_graph.add_bond(idx, jdx)

    return top_graph


# =============================================================================
def find_atomtypes(topology, forcefield, max_iter=10, logger=None):
    """Determine atomtypes for all atoms.

    Parameters
    ----------
    topology : topology.Topology object or TopologyGraph object
    forcefield : foyer.Forcefield
        The forcefield object.
    max_iter : int, optional, default=10
        The maximum number of iterations.
    logger: Log object to write output

    """

    if isinstance(topology, Topology):
        topo_graph = _convert_topology_to_topologygraph(topology)
    elif not isinstance(topology, TopologyGraph):
        m1 = "\t\tERROR: topology must be a Topology or TopologyGraph object in find_atomtypes (atomtyper.py file).\n"
        m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
        m = "\t\t" + len(m1) * "*" + "\n"
        print("\n" + m + m1 + m2 + m) if logger is None else logger.info("\n" + m + m1 + m2 + m)
        exit()
    else:
        topo_graph = topology

    typemap = {
        atom_index: {"whitelist": set(), "blacklist": set(), "atomtype": None}
        for atom_index in topology._graphdict.keys()
    }

    rules = _load_rules(forcefield, typemap)

    # Only consider rules for elements found in topology
    subrules = dict()

    system_elements = set()
    for idx, bonded_list in topology._graphdict.items():
        # First add non-element types, which are strings, then elements
        name = topology._names[idx]
        if name.startswith("_"):
            if name in forcefield.non_element_types:
                system_elements.add(name)
        else:
            element = topology.elements[idx]
            at_number = atomic_number[element]
            at_symbol = element
            try:
                element_from_num = ele.element_from_atomic_number(at_number).symbol
                element_from_sym = ele.element_from_symbol(at_symbol).symbol
                assert element_from_num == element_from_sym
                system_elements.add(element_from_num)
            except ElementError:
                m1 = "\t\tERROR: Atom {} as having neither an element " \
                     "nor non-element type. (atomtyper.py file).\n".format(name)
                m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
                m = "\t\t" + len(m1) * "*" + "\n"
                print("\n" + m + m1 + m2 + m) if logger is None else logger.info("\n" + m + m1 + m2 + m)
                exit()
                pass
            except AssertionError:
                m1 = "\t\tERROR: Atom {} has mismatching atom number {} " \
                     "and atom symbol {}. (atomtyper.py file).\n".format(name, at_number, at_symbol)
                m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
                m = "\t\t" + len(m1) * "*" + "\n"
                print("\n" + m + m1 + m2 + m) if logger is None else logger.info("\n" + m + m1 + m2 + m)
                exit()

    for key, val in rules.items():
        atom = val.nodes[0]["atom"]
        if len(list(atom.find_data("atom_symbol"))) == 1 and not list(atom.find_data("not_expression")):
            try:
                element = next(atom.find_data("atom_symbol")).children[0]
            except IndexError:
                try:
                    atomic_num = next(atom.find_data("atomic_num")).children[0]
                    element = ele.element_from_atomic_number(atomic_num).symbol
                except IndexError:
                    element = None
        else:
            element = None
        if element is None or element in system_elements:
            subrules[key] = val
    rules = subrules

    _iterate_rules(rules, topo_graph, typemap, max_iter=max_iter, logger=None)
    _resolve_atomtypes(topo_graph, typemap, logger=logger)

    return typemap


# =============================================================================
def _load_rules(forcefield, typemap):
    """Load atomtyping rules from a forcefield into SMARTSGraphs."""
    rules = dict()
    # For every SMARTS string in the force field,
    # create a SMARTSGraph object
    for rule_name, smarts in forcefield.atomTypeDefinitions.items():
        if not smarts:  # We want to skip over empty smarts definitions
            continue
        overrides = forcefield.atomTypeOverrides.get(rule_name)
        if overrides is not None:
            overrides = set(overrides)
        else:
            overrides = set()
        rules[rule_name] = SMARTSGraph(
            smarts_string=smarts,
            parser=forcefield.parser,
            name=rule_name,
            overrides=overrides,
            typemap=typemap,
        )
    return rules


# =============================================================================
def _iterate_rules(rules, topology_graph, typemap, max_iter, logger=None):
    """Iterate through all the rules until the white- and blacklists converge.

    Parameters
    ----------
    rules : dict
        A dictionary mapping rule names (typically atomtype names) to
        SMARTSGraphs that evaluate those rules.
    topology_graph : TopologyGraph
        The topology graph that we are trying to atomtype.
    max_iter : int
        The maximum number of iterations.

    """
    for _ in range(max_iter):
        max_iter -= 1
        found_something = False
        for rule in rules.values():
            for match_index in rule.find_matches(topology_graph, typemap):
                atom = typemap[match_index]
                # This conditional is not strictly necessary, but it prevents
                # redundant set addition on later iterations
                if rule.name not in atom["whitelist"]:
                    atom["whitelist"].add(rule.name)
                    atom["blacklist"] |= rule.overrides
                    found_something = True
        if not found_something:
            break
    else:
        m1 = "\t\tWARNING: Reached maximum iterations. Something probably went wrong."
        m = "\t\t" + len(m1) * "*" + "\n"
        print("\n" + m + m1 + m) if logger is None else logger.warning("\n" + m + m1 + m)

    return typemap


def _resolve_atomtypes(topology_graph, typemap, logger=None):
    """Determine the final atomtypes from the white- and blacklists."""
    atoms = {
        atom_idx: data for atom_idx, data in topology_graph.atoms(data=True)
    }
    for atom_id, atom in typemap.items():
        atomtype = [
            rule_name for rule_name in atom["whitelist"] - atom["blacklist"]
        ]
        if len(atomtype) == 1:
            atom["atomtype"] = atomtype[0]
        elif len(atomtype) > 1:
            try:
                m1 = "\t\tERROR: Found multiple types for atom {} ({}): {}.\n".format(
                    atom_id, atoms[atom_id].atomic_number, atomtype)
            except TypeError:
                m1 = "\t\tERROR: Found multiple types for atom {} ({}): {}.\n".format(
                    atom_id, 0, atomtype)

            m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
            m = "\t\t" + len(m1) * "*" + "\n"

            m3 = "\n\t\tAtoms already typed:\n"
            m3 += "\t\t#Index Atomic_Number Type\n"
            for i in range(atom_id):
                try:
                    m3 += "\t\t{0:6d} {1:6d} {2:20s}\n".format(i, atoms[i].atomic_number, typemap[i]['atomtype'])
                except TypeError:
                    m3 += "\t\t{0:6d} {1:6d} {2:20s}\n".format(i, 0, typemap[i]['atomtype'])

            print("\n"+m+m1+m2+m3+m) if logger is None else logger.info("\n"+m+m1+m2+m3+m)
            exit()
        else:
            m1 = "\t\tERROR: Found no types for atom {} ({}).\n".format(atom_id, atoms[atom_id].atomic_number)
            m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
            m = "\t\t" + len(m1) * "*" + "\n"
            m3 = "\n\t\tAtoms already typed:\n"
            m3 += "\t\t#Index Atomic_Number Type\n"
            for i in range(atom_id):
                try:
                    m3 += "\t\t{0:6d} {1:6d} {2:20s}\n".format(i, atoms[i].atomic_number, typemap[i]['atomtype'])
                except TypeError:
                    m3 += "\t\t{0:6d} {1:6d} {2:20s}\n".format(i, 0, typemap[i]['atomtype'])

            print("\n"+m+m1+m2+m3+m) if logger is None else logger.info("\n"+m+m1+m2+m3+m)
            exit()

    m4 = "\n\t\tAll atoms typed:\n"
    m5 = "\t\t" + len(m4) * "*" + "\n"
    m4 += "\t\t#Index Atomic_Number Type\n"

    for i in range(len(atoms)):
        try:
            m4 += "\t\t{0:6d} {1:6d} {2:20s}\n".format(i, atoms[i].atomic_number, typemap[i]['atomtype'])
        except TypeError:
            m4 += "\t\t{0:6d} {1:6d} {2:20s}\n".format(i, 0, typemap[i]['atomtype'])

    print("\n"+m5+m4+m5) if logger is None else logger.info("\n"+m5+m4+m5)
