import os
import re
import glob
import warnings
import topology
import itertools
import numpy as np
import parmed as pmd
import xml.etree.ElementTree as ET
from tempfile import NamedTemporaryFile
from typing import Callable, Iterable, List
from pkg_resources import iter_entry_points, resource_filename
from replicate_polymer_lib.forcefield.Validator import Validator
from replicate_polymer_lib.forcefield.Smart import SMARTS
from replicate_polymer_lib.forcefield import element as elem
from replicate_polymer_lib.forcefield.atomtyper import find_atomtypes
from replicate_polymer_lib.forcefield.topology_graph import TopologyGraph
from topology import Topology
from topology.atomic_data import atomic_number
import openmm as mm
import openmm.unit as u


force_for = {
    mm.app.forcefield.NonbondedGenerator: "NonBondedForce",
    mm.app.forcefield.HarmonicBondGenerator: "HarmonicBondForce",
    mm.app.forcefield.HarmonicAngleGenerator: "HarmonicAngleForce",
    mm.app.forcefield.PeriodicTorsionGenerator: "PeriodicTorsionForce",
    mm.app.forcefield.RBTorsionGenerator: "RBTorsionForce",
}


# =============================================================================================
def preprocess_forcefield_files(forcefield_files=None):
    """Pre-process foyer Forcefield XML files."""
    if forcefield_files is None:
        return None

    preprocessed_files = []

    for xml_file in forcefield_files:
        if not hasattr(xml_file, "read"):
            f = open(xml_file)
            _, suffix = os.path.split(xml_file)
        else:
            f = xml_file
            suffix = ""

        # read and preprocess
        xml_contents = f.read()
        f.close()
        xml_contents = re.sub(
            r"(def\w*=\w*[\"\'])(.*)([\"\'])",
            lambda m: m.group(1)
            + re.sub(r"&(?!amp;)", r"&amp;", m.group(2))
            + m.group(3),
            xml_contents,
        )

        try:
            """
            Sort topology objects by precedence, defined by the number of
            `type` attributes specified, where a `type` attribute indicates
            increased specificity as opposed to use of `class`
            """
            root = ET.fromstring(xml_contents)
            for element in root:
                if "Force" in element.tag:
                    element[:] = sorted(
                        element,
                        key=lambda child: (
                            -1
                            * len(
                                [
                                    attr_name
                                    for attr_name in child.keys()
                                    if "type" in attr_name
                                ]
                            )
                        ),
                    )
            xml_contents = ET.tostring(root, method="xml").decode()
        except ET.ParseError:
            """
            Provide the user with a warning if sorting could not be performed.
            This indicates a bad XML file, which will be passed on to the
            Validator to yield a more descriptive error message.
            """
            warnings.warn(
                "Invalid XML detected. Could not auto-sort topology "
                "objects by precedence."
            )

        # write to temp file
        temp_file = NamedTemporaryFile(suffix=suffix, delete=False)
        with open(temp_file.name, "w") as temp_f:
            temp_f.write(xml_contents)

        # append temp file name to list
        preprocessed_files.append(temp_file.name)

    return preprocessed_files


# =============================================================================================
def get_available_forcefield_loaders() -> List[Callable]:
    """Get a list of available force field loader functions."""
    available_ff_paths = []
    for entry_point in iter_entry_points(group="foyer.forcefields"):
        available_ff_paths.append(entry_point.load())

    return available_ff_paths


# =============================================================================================
def _topology_from_readpdbformat(structure, non_element_types):
    """Convert a ReadPDBFormat Structure to an OpenMM Topology."""
    topology = mm.app.Topology()
    residues = dict()
    for mdanalysis_residue in structure._universe.residues:
        chain = topology.addChain()
        omm_residue = topology.addResidue(mdanalysis_residue.resname, chain)
        # Index residues on name & number, no other info i.e. chain
        residues[(mdanalysis_residue.resname, mdanalysis_residue.resindex)] = omm_residue
    atoms = dict()

    for mdanalysis_atom in structure._universe.atoms:
        name = mdanalysis_atom.name
        if mdanalysis_atom.name in non_element_types:
            element = non_element_types[mdanalysis_atom.name]
        else:
            at_number = atomic_number[mdanalysis_atom.name]
            if isinstance(at_number, int) and at_number != 0:
                element = elem.Element.getByAtomicNumber(at_number)
            else:
                element = elem.Element.getBySymbol(name)

        omm_atom = topology.addAtom(name, element, residues[(mdanalysis_residue.resname, mdanalysis_residue.resindex)])
        omm_atom.id = mdanalysis_atom.type
        atoms[mdanalysis_atom] = omm_atom
        omm_atom.bond_partners = []

    for bond in structure._universe.bonds:
        atom1 = atoms[bond.atoms[0]]
        atom2 = atoms[bond.atoms[1]]
        topology.addBond(atom1, atom2)
        atom1.bond_partners.append(atom2)
        atom2.bond_partners.append(atom1)
    if structure._unitcell is not None:
        topology.setPeriodicBoxVectors(structure._unitcell * mm.unit.angstroms)

    positions = None
    lcords = list()
    for ikey, coord in structure._atom3d_xyz.items():
        #print(coord[0], coord[1], coord[2], mm.Vec3(coord[0], coord[1], coord[2]) * mm.unit.angstroms )
        lcords.append(mm.Vec3(coord[0], coord[1], coord[2]))
    positions = lcords * mm.unit.angstroms

    return topology, positions


# =============================================================================================
def _topology_from_residue(residue, parent_structure):

    """Convert a MDAnalysis Residue to an equivalent Topology_Graph object."""

    topo = Topology()

    for iatom in residue.atoms:
        idx = iatom.index
        topo.add_vertex(idx)
        topo.elements.append('')
        topo._types.append('')

    for iatom in residue.atoms:
        idx = iatom.index
        topo.elements[idx] = iatom.element
        topo._types[idx] = iatom.name
        for jatom in iatom.bonded_atoms:
            jdx = jatom.index
            topo.add_edge((idx, jdx))

    topology_graph = TopologyGraph.from_topomdanalysis_topologygraph(residue)

    return topo, topology_graph


# =============================================================================================
def _check_independent_residues(molecule):

    """Check to see if residues will constitute independent graphs."""
    for ires in molecule._universe.residues:
        atoms_in_residue = set([*ires.atoms])
        bond_partners_in_residue = [
            item
            for sublist in [atom.bonded_atoms for atom in ires.atoms]
            for item in sublist
        ]
        # Handle the case of a 'residue' with no neighbors
        if not bond_partners_in_residue:
            continue
        if set(atoms_in_residue) != set(bond_partners_in_residue):
            return False
        return True


# =============================================================================================
def _unwrap_typemap(structure, residue_map):

    master_typemap = {
        atom.index: {"whitelist": set(), "blacklist": set(), "atomtype": None}
        for atom in structure._universe.atoms
    }
    for res in structure._universe.residues:
        for res_ref, val in residue_map.items():
            # cJ if id(res.name) == id(res_ref):
            if res.resname == res_ref:
                for i, atom in enumerate(res.atoms):
                    master_typemap[int(atom.index)]["atomtype"] = val[i][
                        "atomtype"
                    ]
    return master_typemap


# =============================================================================================
def _separate_urey_bradleys(system, topology):

    """Separate urey bradley bonds from harmonic bonds in OpenMM System.

    Parameters
    ----------
    topology : openmm.app.Topology
        Molecular structure to find atom types of
    system : openmm System

    """
    atoms = [a for a in topology.atoms()]
    bonds = [b for b in topology.bonds()]
    ub_force = mm.HarmonicBondForce()
    harmonic_bond_force = mm.HarmonicBondForce()
    for force_idx, force in enumerate(system.getForces()):
        if isinstance(force, mm.HarmonicBondForce):
            for bond_idx in range(force.getNumBonds()):
                if (
                    atoms[force.getBondParameters(bond_idx)[0]],
                    atoms[force.getBondParameters(bond_idx)[1]],
                ) not in bonds and (
                    atoms[force.getBondParameters(bond_idx)[1]],
                    atoms[force.getBondParameters(bond_idx)[0]],
                ) not in bonds:
                    ub_force.addBond(*force.getBondParameters(bond_idx))
                else:
                    harmonic_bond_force.addBond(
                        *force.getBondParameters(bond_idx)
                    )
            system.removeForce(force_idx)

    system.addForce(harmonic_bond_force)
    system.addForce(ub_force)


# =============================================================================================
def _check_bonds(data, structure, assert_bond_params, logger=None):
    """Check if any bonds lack paramters."""
    if data.bonds:
        missing = [b for b in structure.bonds if b.type is None]
        #cJ
        for ibond in missing:
            print(ibond, ibond.atom1.atom_type, ibond.atom2.atom_type)
        #cJ

        if missing:
            nmissing = len(structure.bonds) - len(missing)
            m1 = "\t\tERROR: Parameters have not been assigned to all bonds.\n " \
                 "\t\tTotal system bonds: {}, Parametrized bonds: {}\n".format(len(structure.bonds), nmissing)
            m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
            m = "\t\t" + len(m1) * "*" + "\n"
            print("\n"+m+m1+m2+m) if logger  is None else logger .info("\n"+m+m1+m2+m)
            exit()


# =============================================================================================
def _check_angles(data, structure, verbose, assert_angle_params, logger=None):
    """Check if all angles were found and parametrized."""
    if verbose:
        for omm_ids in data.angles:
            missing_angle = True
            for pmd_angle in structure.angles:
                pmd_ids = (
                    pmd_angle.atom1.idx,
                    pmd_angle.atom2.idx,
                    pmd_angle.atom3.idx,
                )
                if pmd_ids == omm_ids:
                    missing_angle = False
            if missing_angle:
                m1 = "\t\tERROR: Missing angle with ids {} and types {}.".format(
                        omm_ids, [structure.atoms[idx].type for idx in omm_ids])
                print(m1) if logger is None else logger.info(m1)

        if missing_angle:
            m2 = "\t\tERROR: Molecule cannot be typed!!!!!. Exiting....\n"
            m = "\t\t" + len(m2) * "*" + "\n"
            print("\n" + m + m2 + m) if logger is None else logger.info("\n" + m + m2 + m)
            exit()

    if data.angles and (len(data.angles) != len(structure.angles)):
        m1 = "\t\tERROR: Parameters have not been assigned to all angles.\n" \
             "\t\tTotal system angles: {}, Parameterized angles: {}\n".format(len(data.angles), len(structure.angles))
        m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
        m = "\t\t" + len(m1) * "*" + "\n"
        print("\n" + m + m1 + m2 + m) if logger is None else logger.info("\n" + m + m1 + m2 + m)
        exit()


# =============================================================================================
def _check_dihedrals(data, structure, verbose, assert_dihedral_params, assert_improper_params, logger=None):

    """Check if all dihedrals, including impropers, were found and parametrized."""
    proper_dihedrals = [
        dihedral for dihedral in structure.dihedrals if not dihedral.improper
    ]

    if verbose:
        for omm_ids in data.propers:
            missing_dihedral = True
            for pmd_proper in proper_dihedrals:
                pmd_ids = (
                    pmd_proper.atom1.idx,
                    pmd_proper.atom2.idx,
                    pmd_proper.atom3.idx,
                    pmd_proper.atom4.idx,
                )
                if pmd_ids == omm_ids:
                    missing_dihedral = False
            for pmd_proper in structure.rb_torsions:
                pmd_ids = (
                    pmd_proper.atom1.idx,
                    pmd_proper.atom2.idx,
                    pmd_proper.atom3.idx,
                    pmd_proper.atom4.idx,
                )
                if pmd_ids == omm_ids:
                    missing_dihedral = False
            if missing_dihedral:
                m1 = "\t\tERROR: Missing dihedral with ids {} and types {}.".format(
                        omm_ids, [structure.atoms[idx].type for idx in omm_ids])
                print(m1) if logger is None else logger.info(m1)
        if missing_dihedral:
            m2 = "\t\tERROR: Molecule cannot be typed!!!!!. Exiting....\n"
            m = "\t\t" + len(m2) * "*" + "\n"
            print("\n" + m + m2 + m) if logger is None else logger.info("\n" + m + m2 + m)
            exit()

    if data.propers and len(data.propers) != len(proper_dihedrals) + len(structure.rb_torsions):
        if data.propers and len(data.propers) < len(proper_dihedrals) + len(structure.rb_torsions):
            m1 = "\t\tERROR: Parameters have been assigned to all proper dihedrals.\n" \
                 "\t\tHowever, there are more parameterized dihedrals ({}).\n" \
                 "\t\tthan total system dihedrals ({})." \
                 " This may be due to having multiple periodic dihedrals  for a single dihedral.\n".format\
                (len(proper_dihedrals) + len(structure.rb_torsions), len(data.propers))
            m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
            ll = int(len(m1) / 2)
            m = "\t\t" + ll * "*" + "\n"
            print("\n" + m + m1 + m2 + m) if logger is None else logger.info("\n" + m + m1 + m2 + m)
            exit()
        else:
            m1 = "\t\tERROR: Parameters have not been assigned to all proper dihedrals.\n" \
                 "\t\tTotal system dihedrals: {}, Parameterized dihedrals: {}.\n" \
                 "\t\tNote that if your system contains torsions of Ryckaert-Bellemans functional form," \
                 " all of these torsions are processed as propers.\n".format\
                (len(data.propers),len(proper_dihedrals) + len(structure.rb_torsions))
            m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
            ll = int(len(m1)/2)
            m = "\t\t" + ll * "*" + "\n"
            print("\n" + m + m1 + m2 + m) if logger is None else logger.info("\n" + m + m1 + m2 + m)
            exit()

    improper_dihedrals = [
        dihedral for dihedral in structure.dihedrals if dihedral.improper
    ]
    if data.impropers and len(data.impropers) != len(improper_dihedrals) + len(
        structure.impropers
    ):
        m1 = (
            "\t\tWARN: Parameters have not been assigned to all impropers.\n" \
            "\t\tTotal system impropers: {}, Parameterized impropers: {}.\n " \
            "\t\tNote that if your system contains torsions of Ryckaert-Bellemans functional form,\n"\
            "\t\tall of these torsions are processed as propers.\n".format
            (len(data.impropers), len(improper_dihedrals) + len(structure.impropers))
        )
        m2 = "\t\tThis warn can be ignored if you are sure that impropers are not needed!!!!!\n"
        ll = int(len(m1) / 3)
        m = "\t\t" + ll * "-" + "\n"
        print(m + m1 + m2 + m) if logger is None else logger.info(m + m1 + m2 + m)


# =============================================================================================
class Forcefield(mm.app.ForceField):
    """Forcefield allowing SMARTS based atomtyping inspired by foyer library.
    (https://github.com/mosdef-hub/foyer)

    Parameters
    ----------
    forcefield_files : list of str, optional, default=None
        List of forcefield files to load.


    """
    # =============================================================================================
    def __init__(self, forcefield_files=None, name=None, validation=True, debug=False, logger=None):

        self.atomTypeDefinitions = dict()
        self.atomTypeOverrides = dict()
        self.atomTypeDesc = dict()
        self.atomTypeRefs = dict()
        self.atomTypeClasses = dict()
        self.atomTypeElements = dict()
        self._included_forcefields = dict()
        self.non_element_types = dict()
        self._version = None
        self._name = None
        self._combining_rule = None
        self._logger = logger

        all_files_to_load = []
        if forcefield_files is not None:
            if isinstance(forcefield_files, (list, tuple, set)):
                for file in forcefield_files:
                    all_files_to_load.append(file)
            else:
                all_files_to_load.append(forcefield_files)

        if name is not None:
            try:
                file = self.included_forcefields[name]
            except KeyError:
                raise IOError("Forcefield {} cannot be found".format(name))
            else:
                all_files_to_load.append(file)

        preprocessed_files = preprocess_forcefield_files(all_files_to_load)
        if validation:
            for ff_file_name in preprocessed_files:
                Validator(ff_file_name, debug, logger=logger)

        super(Forcefield, self).__init__(*preprocessed_files)

        if len(preprocessed_files) == 1:
            self._version = self._parse_version_number(preprocessed_files[0])
            self._name = self._parse_name(preprocessed_files[0])
            self._combining_rule = self._parse_combining_rule(preprocessed_files[0])
        elif len(preprocessed_files) > 1:
            self._version = [self._parse_version_number(f) for f in preprocessed_files]
            self._name = [self._parse_name(f) for f in preprocessed_files]
            self._combining_rule = [self._parse_combining_rule(f) for f in preprocessed_files]
            if len(set(self._combining_rule)) == 1:
                self._combining_rule = self._combining_rule[0]
            else:
                m1 = "\t\tERROR: Inconsistent combining_rule among loaded forcefield files\n"
                m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
                m = "\t\t" + len(m1) * "*" + "\n"
                print("\n" + m + m1 + m2 + m) if self._logger is None else self._logger.info("\n" + m + m1 + m2 + m)
                exit()

        for fp in preprocessed_files:
            os.remove(fp)

        self.parser = SMARTS(self.non_element_types, logger=logger)
        self._system_data = None

    # =============================================================================================
    @property
    def version(self):
        """Return version number of the force field XML file."""
        return self._version

    # =============================================================================================
    @property
    def name(self):
        """Return the name of the force field XML."""
        return self._name

    # =============================================================================================
    @property
    def combining_rule(self):
        """Return the combining rule of this force field."""
        return self._combining_rule

    # =============================================================================================
    @property
    def included_forcefields(self):
        """Return a dictionary mappying for all included forcefields."""
        if any(self._included_forcefields):
            return self._included_forcefields

        ff_dir = resource_filename("foyer", "forcefields")
        ff_filepaths = set(glob.glob(os.path.join(ff_dir, "xml/*.xml")))

        for ff_filepath in ff_filepaths:
            _, ff_file = os.path.split(ff_filepath)
            basename, _ = os.path.splitext(ff_file)
            self._included_forcefields[basename] = ff_filepath
        return self._included_forcefields

    # =============================================================================================
    @property
    def lj14scale(self):
        """Get LJ 1-4 scale for this forcefield."""
        try:
            non_bonded_force_gen = self.get_generator(ff=self, gen_type=mm.app.forcefield.NonbondedGenerator)
        except Exception:
            raise AttributeError(
                "Cannot get lj14Scale for the forcefield "
                "because it doesn't have NonBondedForce."
            )
        return non_bonded_force_gen.lj14scale

    # =============================================================================================
    @property
    def coulomb14scale(self):
        """Get Coulombic 1-4 scale for this forcefield."""
        try:
            non_bonded_force_gen = self.get_generator(
                ff=self, gen_type=mm.app.forcefield.NonbondedGenerator
            )
        except Exception:
            raise AttributeError(
                "Cannot get coulomb14scale for the Forcefield "
                "because it doesn't have NonBondedForce."
            )
        return non_bonded_force_gen.coulomb14scale

    # =============================================================================================
    def _parse_version_number(self, forcefield_file):
        with open(forcefield_file, "r") as f:
            tree = ET.parse(f)
            root = tree.getroot()
            try:
                return root.attrib["version"]
            except KeyError:
                m1 = "\t\tWARN: No force field version number found in force field XML file.\n"
                m = "\t\t" + len(m1) * "-" + "\n"
                print(m + m1 + m) if self._logger is None else self._logger.info(m + m1 + m)
                return None

    # =============================================================================================
    def _parse_name(self, forcefield_file):
        with open(forcefield_file, "r") as f:
            tree = ET.parse(f)
            root = tree.getroot()
            try:
                return root.attrib["name"]
            except KeyError:
                m1 = "\t\tWARN: No force field name found in force field XML file.\n"
                m = "\t\t" + len(m1) * "-" + "\n"
                print(m + m1 + m) if self._logger is None else self._logger.info(m + m1 + m)
                return None

    # =============================================================================================
    def _parse_combining_rule(self, forcefield_file):
        with open(forcefield_file, "r") as f:
            tree = ET.parse(f)
            root = tree.getroot()
            try:
                return root.attrib["combining_rule"]
            except KeyError:
                m1 = "\t\tWARN: No combining rule found in force field XML file. Using Lorentz by default.\n"
                m = "\t\t" + len(m1) * "-" + "\n"
                print(m + m1 + m) if self._logger is None else self._logger.info(m + m1 + m)
                return "lorentz"

    # =============================================================================================
    def _create_element(self, element, mass):
        if not isinstance(element, elem.Element):
            try:
                element = elem.get_by_symbol(element)
            except KeyError:
                # Enables support for non-atomistic "element types"
                if element not in self.non_element_types:
                    warnings.warn(
                        "Non-atomistic element type detected. "
                        "Creating custom element for {}".format(element)
                    )
                element = elem.Element(
                    number=0, mass=mass, name=element, symbol=element
                )
            else:
                return element, False

        return element, True

    # =============================================================================================
    def registerAtomType(self, parameters):
        """Register a new atom type."""
        name = parameters["name"]
        if name in self._atomTypes:
            raise ValueError(
                "Found multiple definitions for atom type: " + name
            )
        atom_class = parameters["class"]
        mass = mm.app.forcefield._convertParameterToNumber(parameters["mass"])
        element = None
        if "element" in parameters:
            element, custom = self._create_element(parameters["element"], mass)
            if custom:
                self.non_element_types[element.symbol] = element

        self._atomTypes[name] = self.__class__._AtomType(
            name, atom_class, mass, element
        )
        if atom_class in self._atomClasses:
            type_set = self._atomClasses[atom_class]
        else:
            type_set = set()
            self._atomClasses[atom_class] = type_set
        type_set.add(name)
        self._atomClasses[""].add(name)

        name = parameters["name"]
        if "def" in parameters:
            self.atomTypeDefinitions[name] = parameters["def"]
        if "overrides" in parameters:
            overrides = set(
                atype.strip() for atype in parameters["overrides"].split(",")
            )
            if overrides:
                self.atomTypeOverrides[name] = overrides
        if "desc" in parameters:
            self.atomTypeDesc[name] = parameters["desc"]
        if "doi" in parameters:
            dois = set(doi.strip() for doi in parameters["doi"].split(","))
            self.atomTypeRefs[name] = dois
        if "element" in parameters:
            self.atomTypeElements[name] = parameters["element"]
        if "class" in parameters:
            self.atomTypeClasses[name] = parameters["class"]

    # =============================================================================================
    def apply(self, molecule, references_file=None, use_residue_map=True,
              assert_bond_params=True, assert_angle_params=True,
              assert_dihedral_params=True, assert_improper_params=False,
              verbose=False, *args, **kwargs):

        """Apply the force field to a molecular structure.

        Parameters
        ----------
        structure : parmed.Structure or mbuild.Compound
            Molecular structure to apply the force field to
        references_file : str, optional, default=None
            Name of file where force field references will be written (in Bibtex
            format)
        use_residue_map : boolean, optional, default=True
            Store atomtyped topologies of residues to a dictionary that maps
            them to residue names.  Each topology, including atomtypes, will be
            copied to other residues with the same name. This avoids repeatedly
            calling the subgraph isomorphism on idential residues and should
            result in better performance for systems with many identical
            residues, i.e. a box of water. Note that for this to be applied to
            independent molecules, they must each be saved as different
            residues in the topology.
        assert_bond_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all system
            bonds.
        assert_angle_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all system
            angles.
        assert_dihedral_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all system
            proper dihedrals.
        assert_improper_params : bool, optional, default=False
            If True, Foyer will exit if parameters are not found for all system
            improper dihedrals.
        verbose : bool, optional, default=False
            If True, Foyer will print debug-level information about notable or
            potentially problematic details it encounters.
        """

        if self.atomTypeDefinitions == {}:
            m1 = "\t\tERROR: Attempting to atom-type using a force field with no atom type defitions.\n"
            m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
            m = "\t\t" + len(m1) * "*" + "\n"
            print("\n"+m+m1+m2+m) if self._logger is None else self._logger.info("\n"+m+m1+m2+m)
            exit()

        if not isinstance(molecule, topology.readmol.readpdbformat.ReadPdbFormat):
            m1 = "\t\tERROR: Molecule must be a ReadPdbFormat object. Considerer to use a PDB file as input.\n"
            m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
            m = "\t\t" + len(m1) * "*" + "\n"
            print("\n"+m+m1+m2+m) if self._logger is None else self._logger.info("\n"+m+m1+m2+m)
            exit()

        typemap = self.run_atomtyping(molecule, use_residue_map=use_residue_map, **kwargs)

        self._apply_typemap(molecule, typemap)

        pmd_structure = self.parametrize_system(structure=molecule, references_file=references_file,
                                                assert_bond_params=assert_bond_params, assert_angle_params=assert_angle_params,
                                                assert_dihedral_params=assert_dihedral_params,
                                                assert_improper_params=assert_improper_params,
                                                verbose=verbose, *args, **kwargs,)

        return pmd_structure

    # =============================================================================================
    def _create_element(self, element, mass):
        if not isinstance(element, elem.Element):
            try:
                element = elem.get_by_symbol(element)
            except KeyError:
                # Enables support for non-atomistic "element types"
                if element not in self.non_element_types:
                    warnings.warn(
                        "Non-atomistic element type detected. "
                        "Creating custom element for {}".format(element)
                    )
                element = elem.Element(
                    number=0, mass=mass, name=element, symbol=element
                )
            else:
                return element, False

        return element, True

    # =============================================================================================
    def run_atomtyping(self, structure, use_residue_map=True, **kwargs):
        """Atomtype the topology.

        Parameters
        ----------
        structure : Molecular structure ReadPdbFormat
            Molecular structure to find atom types of
        use_residue_map : boolean, optional, default=True
            Store atomtyped topologies of residues to a dictionary that maps
            them to residue names.  Each topology, including atomtypes, will be
            copied to other residues with the same name. This avoids repeatedly
            calling the subgraph isomorphism on idential residues and should
            result in better performance for systems with many identical
            residues, i.e. a box of water. Note that for this to be applied to
            independent molecules, they must each be saved as different
            residues in the topology.
        """

        if use_residue_map:
            independent_residues = _check_independent_residues(structure)
            if independent_residues:
                residue_map = dict()
                # Need to call this only once and store results for later id() comparisons
                for res_id, res in enumerate(structure._universe.residues):
                    if (
                        structure._universe.residues[res_id].resname
                        not in residue_map.keys()
                    ):
                        topo_res, topograph_res = _topology_from_residue(res, structure)
                        typemap = find_atomtypes(topo_res, forcefield=self)
                        residue_map[res.resname] = typemap

                typemap = _unwrap_typemap(structure, residue_map)

            else:

                typemap = find_atomtypes(structure._topology, forcefield=self)

        else:
            typemap = find_atomtypes(structure._topology, forcefield=self)

        return typemap

    # =============================================================================================
    def parametrize_system(self, structure=None, references_file=None, assert_bond_params=True,
                           assert_angle_params=True, assert_dihedral_params=True,
                           assert_improper_params=False, verbose=False, *args, **kwargs):

        """Create system based on resulting typemapping"""
        topology, positions = _topology_from_readpdbformat(structure, self.non_element_types)

        system = self.createSystem(topology, *args, **kwargs)

        _separate_urey_bradleys(system, topology)

        data = self._system_data

        structure = pmd.openmm.load_topology(topology=topology, system=system)
        structure.bonds.sort(key=lambda x: x.atom1.idx)
        structure.positions = positions
        box_vectors = topology.getPeriodicBoxVectors()
        if box_vectors is not None:
            structure.box_vectors = box_vectors

        _check_bonds(data, structure, assert_bond_params, logger=self._logger)
        _check_angles(data, structure, verbose, assert_angle_params, logger=self._logger)
        _check_dihedrals(data, structure, verbose, assert_dihedral_params, assert_improper_params, logger=self._logger)

        if references_file:
            atom_types = set(atom.type for atom in structure.atoms)
            self._write_references_to_file(atom_types, references_file)

        try:
            structure.combining_rule = self.combining_rule
        except ValueError as e:
            m1 = "\t\tERROR: Combination rule {} is not implemented must " \
                 "be 'lorentz' or 'geometric'.\n".format(self.combining_rule)
            m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
            m = "\t\t" + len(m1) * "*" + "\n"
            print("\n"+m+m1+m2+m) if self._logger is None else self._logger.info("\n"+m+m1+m2+m)
            exit()

        if self.combining_rule == "geometric":
            self._patch_parmed_adjusts(structure, combining_rule=self.combining_rule)

        total_charge = sum([atom.charge for atom in structure.atoms])
        if not np.allclose(total_charge, 0.0):
            m1 = ("\t\tWARN: Parametrized structure has non-zero charge.\n" \
                  "\t\tStructure's total charge: {}".format(total_charge))
            m2 = "\t\tCheck charges!!!!!\n"
            ll = int(len(m1) / 1)
            m = "\t\t" + ll * "-" + "\n"
            print(m + m1 + m2 + m) if self._logger is None else self._logger.info(m + m1 + m2 + m)

        return structure

    # =============================================================================================
    def createSystem(self, topology, nonbondedMethod=mm.app.forcefield.NoCutoff,
                     nonbondedCutoff=1.0 * mm.unit.nanometer, constraints=None, rigidWater=True,
                     removeCMMotion=True, hydrogenMass=None, switchDistance=None, **args):

        """Construct an OpenMM System representing a Topology with this force field.

        Parameters
        ----------
        topology : Topology
            The Topology for which to create a System
        nonbondedMethod : object=NoCutoff
            The method to use for nonbonded interactions.  Allowed values are
            NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, or PME.
        nonbondedCutoff : distance=1*nanometer
            The cutoff distance to use for nonbonded interactions
        constraints : object=None
            Specifies which bonds and angles should be implemented with constraints.
            Allowed values are None, HBonds, AllBonds, or HAngles.
        rigidWater : boolean=True
            If true, water molecules will be fully rigid regardless of the value
            passed for the constraints argument
        removeCMMotion : boolean=True
            If true, a CMMotionRemover will be added to the System
        hydrogenMass : mass=None
            The mass to use for hydrogen atoms bound to heavy atoms.  Any mass
            added to a hydrogen is subtracted from the heavy atom to keep
            their total mass the same.
        switchDistance : float=None
            The distance at which the potential energy switching function is turned on for
        args
             Arbitrary additional keyword arguments may also be specified.
             This allows extra parameters to be specified that are specific to
             particular force fields.

        Returns
        -------
        system
            the newly created System
        """
        args["switchDistance"] = switchDistance
        # Overwrite previous _SystemData object
        data = mm.app.ForceField._SystemData(topology)
        self._system_data = data

        # TODO: Better way to lookup nonbonded parameters...?
        nonbonded_params = None
        for generator in self.getGenerators():
            if isinstance(generator, mm.app.forcefield.NonbondedGenerator):
                nonbonded_params = generator.params.paramsForType
                break

        for chain in topology.chains():
            for res in chain.residues():
                for atom in res.atoms():
                    data.atomType[atom] = atom.id
                    if nonbonded_params:
                        params = nonbonded_params[atom.id]
                        data.atomParameters[atom] = params

        # Create the System and add atoms
        sys = mm.System()
        for atom in topology.atoms():
            # Look up the atom type name, returning a helpful error message if it cannot be found.
            if atom not in data.atomType:
                raise Exception(
                    "Could not identify atom type for atom '%s'." % str(atom)
                )
            typename = data.atomType[atom]

            # Look up the type name in the list of registered atom types,
            # returning a helpful error message if it cannot be found.
            if typename not in self._atomTypes:
                msg = (
                    "Could not find typename '%s' for atom '%s' in list of known atom types.\n"
                    % (typename, str(atom))
                )
                msg += "Known atom types are: %s" % str(self._atomTypes.keys())
                raise Exception(msg)

            # Add the particle to the OpenMM system.
            mass = self._atomTypes[typename].mass
            sys.addParticle(mass)

        # Adjust hydrogen masses if requested.
        if hydrogenMass is not None:
            if not u.is_quantity(hydrogenMass):
                hydrogenMass *= mm.unit.dalton
            for atom1, atom2 in topology.bonds():
                if atom1.element == elem.hydrogen:
                    (atom1, atom2) = (atom2, atom1)
                if atom2.element == elem.hydrogen and atom1.element not in (
                    elem.hydrogen,
                    None,
                ):
                    transfer_mass = hydrogenMass - sys.getParticleMass(
                        atom2.index
                    )
                    sys.setParticleMass(atom2.index, hydrogenMass)
                    mass = sys.getParticleMass(atom1.index) - transfer_mass
                    sys.setParticleMass(atom1.index, mass)

        # Set periodic boundary conditions.
        box_vectors = topology.getPeriodicBoxVectors()
        if box_vectors is not None:
            sys.setDefaultPeriodicBoxVectors(
                box_vectors[0], box_vectors[1], box_vectors[2]
            )
        elif nonbondedMethod not in [mm.app.forcefield.NoCutoff, mm.app.forcefield.CutoffNonPeriodic]:
            raise ValueError(
                "Requested periodic boundary conditions for a "
                "Topology that does not specify periodic box "
                "dimensions"
            )

        # Make a list of all unique angles
        unique_angles = set()
        for bond in data.bonds:
            for atom in data.bondedToAtom[bond.atom1]:
                if atom != bond.atom2:
                    if atom < bond.atom2:
                        unique_angles.add((atom, bond.atom1, bond.atom2))
                    else:
                        unique_angles.add((bond.atom2, bond.atom1, atom))
            for atom in data.bondedToAtom[bond.atom2]:
                if atom != bond.atom1:
                    if atom > bond.atom1:
                        unique_angles.add((bond.atom1, bond.atom2, atom))
                    else:
                        unique_angles.add((atom, bond.atom2, bond.atom1))
        data.angles = sorted(list(unique_angles))

        # Make a list of all unique proper torsions
        unique_propers = set()
        for angle in data.angles:
            for atom in data.bondedToAtom[angle[0]]:
                if atom not in angle:
                    if atom < angle[2]:
                        unique_propers.add((atom, angle[0], angle[1], angle[2]))
                    else:
                        unique_propers.add((angle[2], angle[1], angle[0], atom))
            for atom in data.bondedToAtom[angle[2]]:
                if atom not in angle:
                    if atom > angle[0]:
                        unique_propers.add((angle[0], angle[1], angle[2], atom))
                    else:
                        unique_propers.add((atom, angle[2], angle[1], angle[0]))
        data.propers = sorted(list(unique_propers))

        # Make a list of all unique improper torsions
        for atom in range(len(data.bondedToAtom)):
            bonded_to = data.bondedToAtom[atom]
            if len(bonded_to) > 2:
                for subset in itertools.combinations(bonded_to, 3):
                    data.impropers.append(
                        (atom, subset[0], subset[1], subset[2])
                    )

        # Identify bonds that should be implemented with constraints
        if constraints == mm.app.forcefield.AllBonds or constraints == mm.app.forcefield.HAngles:
            for bond in data.bonds:
                bond.isConstrained = True
        elif constraints == mm.app.forcefield.HBonds:
            for bond in data.bonds:
                atom1 = data.atoms[bond.atom1]
                atom2 = data.atoms[bond.atom2]
                bond.isConstrained = atom1.name.startswith(
                    "H"
                ) or atom2.name.startswith("H")
        if rigidWater:
            for bond in data.bonds:
                atom1 = data.atoms[bond.atom1]
                atom2 = data.atoms[bond.atom2]
                if atom1.residue.name == "HOH" and atom2.residue.name == "HOH":
                    bond.isConstrained = True

        # Identify angles that should be implemented with constraints
        if constraints == mm.app.forcefield.HAngles:
            for angle in data.angles:
                atom1 = data.atoms[angle[0]]
                atom2 = data.atoms[angle[1]]
                atom3 = data.atoms[angle[2]]
                numH = 0
                if atom1.name.startswith("H"):
                    numH += 1
                if atom3.name.startswith("H"):
                    numH += 1
                data.isAngleConstrained.append(
                    numH == 2 or (numH == 1 and atom2.name.startswith("O"))
                )
        else:
            data.isAngleConstrained = len(data.angles) * [False]
        if rigidWater:
            for i in range(len(data.angles)):
                angle = data.angles[i]
                atom1 = data.atoms[angle[0]]
                atom2 = data.atoms[angle[1]]
                atom3 = data.atoms[angle[2]]
                if (
                    atom1.residue.name == "HOH"
                    and atom2.residue.name == "HOH"
                    and atom3.residue.name == "HOH"
                ):
                    data.isAngleConstrained[i] = True

        # Add virtual sites
        for atom in data.virtualSites:
            (site, atoms, excludeWith) = data.virtualSites[atom]
            index = atom.index
            data.excludeAtomWith[excludeWith].append(index)
            if site.type == "average2":
                sys.setVirtualSite(
                    index,
                    mm.TwoParticleAverageSite(
                        atoms[0], atoms[1], site.weights[0], site.weights[1]
                    ),
                )
            elif site.type == "average3":
                sys.setVirtualSite(
                    index,
                    mm.ThreeParticleAverageSite(
                        atoms[0],
                        atoms[1],
                        atoms[2],
                        site.weights[0],
                        site.weights[1],
                        site.weights[2],
                    ),
                )
            elif site.type == "outOfPlane":
                sys.setVirtualSite(
                    index,
                    mm.OutOfPlaneSite(
                        atoms[0],
                        atoms[1],
                        atoms[2],
                        site.weights[0],
                        site.weights[1],
                        site.weights[2],
                    ),
                )
            elif site.type == "localCoords":
                local_coord_site = mm.LocalCoordinatesSite(
                    atoms[0],
                    atoms[1],
                    atoms[2],
                    mm.Vec3(
                        site.originWeights[0],
                        site.originWeights[1],
                        site.originWeights[2],
                    ),
                    mm.Vec3(
                        site.xWeights[0], site.xWeights[1], site.xWeights[2]
                    ),
                    mm.Vec3(
                        site.yWeights[0], site.yWeights[1], site.yWeights[2]
                    ),
                    mm.Vec3(
                        site.localPos[0], site.localPos[1], site.localPos[2]
                    ),
                )
                sys.setVirtualSite(index, local_coord_site)

        # Add forces to the System
        for force in self._forces:
            force.createForce(sys, data, nonbondedMethod, nonbondedCutoff, args)
        if removeCMMotion:
            sys.addForce(mm.CMMotionRemover())

        # Let force generators do postprocessing
        for force in self._forces:
            if "postprocessSystem" in dir(force):
                force.postprocessSystem(sys, data, args)

        # Execute scripts found in the XML files.
        for script in self._scripts:
            exec(script, locals())

        return sys

    # =============================================================================================
    def _apply_typemap(self, structure, typemap):

        """Add atomtypes to the topology according to the typemap."""
        for atom in structure._universe.atoms:
            atom.type = typemap[atom.index]["atomtype"]

        if not all([a.type for a in structure._universe.atoms]):
            m1 = "\t\tERROR: Not all atoms in topology have atom types. (forcefield.py: _apply_typemap).\n"
            m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
            m = "\t\t" + len(m1) * "*" + "\n"
            print("\n"+m+m1+m2+m) if self._logger is None else self._logger.info("\n"+m+m1+m2+m)
            exit()

    # =============================================================================================
    @staticmethod
    def _patch_parmed_adjusts(structure, combining_rule="geometric"):
        if combining_rule != "geometric":
            raise BaseException()

        for adj in structure.adjusts:
            sig1 = adj.atom1.sigma
            sig2 = adj.atom2.sigma
            sig_geo = (sig1 * sig2) ** 0.5
            adj.type.sigma = sig_geo

        return structure













