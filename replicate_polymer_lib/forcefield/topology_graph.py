"""Module to represent chemical systems as graph structures."""
import networkx as nx
from topology.atomic_data import atomic_number


class AtomData:
    """Store atom data necessary for atom typing.

    Parameters
    ----------
    index: int
        The index of the atom in the topology graph
    name: str
        The name of the atom
    atomic_number: optional, int, default=None
        The atomic number if the atom represents an element
    element: optional, str, default=None
        The element symbol associated with the atom
        if it represents an element
    **kwargs
        Any extra information for this atom

    ToDo: After 3.6 support removed, convert this to a dataclass
    """

    def __init__(self, index, name, atomic_number=None, element=None, **kwargs):
        self.index = index
        self.name = name
        self.atomic_number = atomic_number
        self.element = element
        for key, value in kwargs.items():
            setattr(self, key, value)


class TopologyGraph(nx.Graph):
    """A general TopologyGraph.

    This class subclasses nx.Graph to provide a general
    topology Graph for atom typing in foyer. Each node in
    this graph is identified by atom indices and edge in this
    graph represents a bond between the atoms.
    """

    def __init__(self, *args, **kwargs):
        super(TopologyGraph, self).__init__(*args, **kwargs)

    def add_atom(self, index, name, atomic_number=None, element=None, logger=None, **kwargs):
        """Add an atom to the topology graph.

        Parameters
        ----------
        index: int
            The index of the atom in the topology graph
        name: str
            The name of the atom. For non-element type,
            the name must start with an underscore (_)
        atomic_number: optional, int, default=None
            The atomic number if the atom represents an element
        element: optional, str, default=None
            The element symbol associated with the atom
            if it represents an element
        **kwargs
            Any extra information for this atom

        See Also
        --------
        foyer.topology_graph.AtomData
            The class used to store atom data
        """
        if not name.startswith("_") and not (atomic_number and element):
            m1 = "\t\tERROR: For atoms representing an element, please include " \
                 "both the atomic_number or element symbol for the atom\n"
            m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
            m = "\t\t" + len(m1) * "*" + "\n"
            print("\n" + m + m1 + m2 + m) if logger is None else logger.info("\n" + m + m1 + m2 + m)
            exit()

        atom_data = AtomData(index, name, atomic_number, element, **kwargs)
        self.add_node(index, atom_data=atom_data)

    def add_bond(self, atom_1_index, atom_2_index):
        """Add a bond(edge) between two atoms in this TopologyGraph.

        Parameters
        ----------
        atom_1_index: int
            The index of the first atom that forms this bond
        atom_2_index: int
            The index of the second atom that forms this bond
        """
        self.add_edge(atom_1_index, atom_2_index)

    def atoms(self, data=False):
        """Iterate through atoms in the TopologyGraph."""
        if data:
            for idx, data in self.nodes(data=data):
                yield idx, data["atom_data"]
        else:
            for idx in self.nodes(data=data):
                yield idx

    def add_bond_partners(self):
        """Add atom indices for atoms involved in a bond."""
        for atom_idx, data in self.nodes(data=True):
            data["bond_partners"] = list(self.neighbors(atom_idx))

    @classmethod
    def from_gmso_topology(cls, gmso_topology):
        """Return a TopologyGraph with relevant attributes from an GMSO topology.

        Parameters
        ----------
        gmso_topology: gmso.Topology
            The GMSO Topology

        Returns
        -------
        TopologyGraph
            The equivalent TopologyGraph of the openFF Topology `openff_topology`
        """
        from replicate_polymer_lib.forcefield.utils.io import import_

        gmso = import_("gmso")  # This might only be required for type checking

        if not isinstance(gmso_topology, gmso.Topology):
            raise TypeError(
                f"Expected `openff_topology` to be of type {gmso.Topology}. "
                f"Got {type(gmso_topology).__name__} instead"
            )

        top_graph = cls()
        for atom in gmso_topology.sites:
            if isinstance(atom, gmso.Atom):
                if atom.name.startswith("_"):
                    top_graph.add_atom(
                        name=atom.name,
                        index=gmso_topology.get_index(atom),
                        atomic_number=None,
                        element=atom.name,
                    )

                else:
                    top_graph.add_atom(
                        name=atom.name,
                        index=gmso_topology.get_index(atom),
                        atomic_number=atom.element.atomic_number,
                        element=atom.element.symbol,
                    )

        for top_bond in gmso_topology.bonds:
            atoms_indices = [
                gmso_topology.get_index(atom)
                for atom in top_bond.connection_members
            ]
            top_graph.add_bond(atoms_indices[0], atoms_indices[1])

        return top_graph

    @classmethod
    def from_topomdanalysis_topologygraph(cls, residue_mdtraj, logger=None):

        """
        J.Ramos

        Create a topologygraph object from a Residue
        instance of MDAnalysis (https://www.mdanalysis.org/)

        Args:
            residue_mdtraj: A Mdtraj residue

        Returns:
            A topology_graph instance.

        """

        import MDAnalysis

        top_graph = cls()
        if  not isinstance(residue_mdtraj, MDAnalysis.core.groups.Residue):
            m1 = "\t\tERROR: Expected `topology` to be of type `MDAnalysis.core.topology.Residue`. " \
                 "\n\t\tGot {} instead. (topology_graph.py: from_topomdanalysis_topologygraph)\n".format(type(residue_mdtraj))
            m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
            m = "\t\t" + len(m1) * "*" + "\n"
            print("\n" + m + m1 + m2 + m) if logger is None else logger.info("\n" + m + m1 + m2 + m)
            exit()

        for iatom in residue_mdtraj.atoms:
            # For no-real atoms, rather interaction sites (i.e: CG or UA models)
            if iatom.name.startswith("_"):
                top_graph.add_atom(name=iatom.name,
                                   index=iatom.index,
                                   atomic_number=None,
                                   element=iatom.name)
            else:
                top_graph.add_atom(name=iatom.name,
                                   index=iatom.index,
                                   atomic_number=atomic_number[iatom.element],
                                   element=iatom.element)

        for iatom in residue_mdtraj.atoms:
            idx = iatom.index
            for jatom in iatom.bonded_atoms:
                jdx = jatom.index
                if jdx > idx:
                    top_graph.add_bond(jdx, idx)
                else:
                    top_graph.add_bond(idx, jdx)

        return top_graph
