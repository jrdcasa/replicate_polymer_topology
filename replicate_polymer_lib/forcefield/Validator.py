"""Validaton of force field XML files."""
import os
import datetime
import lark
import networkx as nx

from collections import Counter
from os.path import abspath, join, split
from lxml import etree
from lxml.etree import DocumentInvalid

from replicate_polymer_lib.forcefield.Smarts_Graph import SMARTSGraph


class Validator(object):
    """Verifies formatting of force field XML."""

    # =============================================================================================
    def __init__(self, ff_file_name, debug=False, logger=None):

        from replicate_polymer_lib.forcefield.Forcefield import preprocess_forcefield_files

        self._logger = logger

        try:
            preprocessed_ff_file_name = preprocess_forcefield_files([ff_file_name])

            self.ff_tree = etree.parse(preprocessed_ff_file_name[0])
            self.validate_xsd(self.ff_tree)

            self.atom_type_names = self.ff_tree.xpath("/ForceField/AtomTypes/Type/@name")
            self.atom_types = self.ff_tree.xpath("/ForceField/AtomTypes/Type")

            self.validate_class_type_exclusivity(self.ff_tree)

            # Loading forcefield should succeed, because XML can be parsed and
            # basics have been validated.
            from replicate_polymer_lib.forcefield.Forcefield import Forcefield
            self.smarts_parser = Forcefield(preprocessed_ff_file_name, validation=False).parser

        finally:
            for ff_file_name in preprocessed_ff_file_name:
                os.remove(ff_file_name)

        self.validate_smarts(debug=debug)
        self.validate_overrides()
        self.all_forces_included_xml = self.get_all_forces_cj()

    # =============================================================================================
    def validate_xsd(self, ff_tree, xsd_file=None):
        """Check consistency with forcefields/ff.xsd."""
        if xsd_file is None:
            xsd_file = join(split(abspath(__file__))[0], "../../forcefields", "ff.xsd")

        try:
            xmlschema_doc = etree.parse(xsd_file)
        except OSError:
            now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
            m1 = "\t\tERROR: XMLScheme XSD file for XML is missing {} ({})\n".format(xsd_file, now)
            m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
            m = "\t\t" + len(m1) * "*" + "\n"
            print("\n"+m+m1+m2+m) if self._logger is None else self._logger.info("\n"+m+m1+m2+m)
            exit()

        try:
            xmlschema = etree.XMLSchema(xmlschema_doc)
        except etree.XMLSchemaParseError as e:
            now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
            m1 = "\t\tERROR: {}\n".format(e)
            m1 += "\t\tXSD_File = {}\n".format(xsd_file)
            m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
            m = "\t\t" + len(m1) * "*" + "\n"
            print("\n"+m+m1+m2+m) if self._logger is None else self._logger.info("\n"+m+m1+m2+m)
            exit()

        error_texts = {
            "missing_atom_type_in_nonbonded": (
                "Atom type {} is found in NonbondedForce at line {}"
                " but undefined in AtomTypes"
            ),
            "nonunique_atomtype_name": "Atom type {} is defined a second time at line {}",
            "atomtype_name_key": "Atom type {} is defined a second time at line {}",
        }

        try:
            xmlschema.assertValid(ff_tree)
        except DocumentInvalid as ex:
            message = ex.error_log.last_error.message
            line = ex.error_log.last_error.line
            now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
            m1 = "\t\tERROR: {} in line {}\n".format(message, line)
            m1 += "\t\tXSD_File = {}\n".format(xsd_file)
            m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
            m = "\t\t" + len(m1) * "*" + "\n"
            print("\n"+m+m1+m2+m) if self._logger is None else self._logger.info("\n"+m+m1+m2+m)
            exit()

    # =============================================================================================
    def validate_class_type_exclusivity(self, ff_tree):
        """Assert unique bond/angle/dihedral definitions and class lengths."""
        sections = {
            "HarmonicBondForce/Bond": 2,
            "HarmonicAngleForce/Angle": 3,
            "RBTorsionForce/Proper": 4,
        }

        errors = []
        for element, num_atoms in sections.items():
            valid_attribs = set()
            for n in range(1, num_atoms + 1):
                valid_attribs.add("class{}".format(n))
                valid_attribs.add("type{}".format(n))

            for entry in ff_tree.xpath("/ForceField/{}".format(element)):
                attribs = [
                    valid
                    for valid in valid_attribs
                    if entry.attrib.get(valid) is not None
                ]
                if num_atoms != len(attribs):
                    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
                    m1 = '\t\tERROR: Invalid number of "class" and/or "type" attributes for \n' \
                         '\t\t\t{} at line {}\n'.format(element, entry.sourceline)
                    m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
                    m = "\t\t" + len(m1) * "*" + "\n"
                    print("\n" + m + m1 + m2 + m) if self._logger is None else self._logger.info("\n" + m + m1 + m2 + m)
                    exit()

                number_endings = Counter([a[-1] for a in attribs])
                if not all(1 == x for x in number_endings.values()):
                    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
                    m1 = '\t\tERROR: Only one "class" or "type" attribute may be defined for ' \
                         'each atom in a bonded force\n' \
                         '\t\t\tat line {}\n'.format(entry.sourceline)
                    m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
                    m = "\t\t" + len(m1) * "*" + "\n"
                    print("\n" + m + m1 + m2 + m) if self._logger is None else self._logger.info("\n" + m + m1 + m2 + m)
                    exit()

                referenced_types = []
                for valid in valid_attribs:
                    if valid.startswith("type"):
                        atomtype = entry.attrib.get(valid)
                        if atomtype:
                            referenced_types.append(atomtype)

                for atomtype in referenced_types:
                    if atomtype not in self.atom_type_names:
                        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
                        m1 = '\t\tERROR: Atom type {} is found in {} at line {} but ' \
                             ' undefined in AtomTypes\n'.format(atomtype, element.split("/")[0], entry.sourceline)
                        m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
                        m = "\t\t" + len(m1) * "*" + "\n"
                        print("\n" + m + m1 + m2 + m) if self._logger is None else self._logger.info(
                            "\n" + m + m1 + m2 + m)
                        exit()

    # =============================================================================================
    def validate_smarts(self, debug):
        """Check SMARTS definitions for missing or non-parseable."""
        missing_smarts = []
        errors = []
        for entry in self.atom_types:
            smarts_string = entry.attrib.get("def")
            if not smarts_string and debug:
                m = "\t\t WARN: You have empty smart definition(s) for {}".format(entry.attrib['name'])
                print(m) if self._logger is None else self._logger.info(m)
                continue
            name = entry.attrib["name"]
            if smarts_string is None:
                missing_smarts.append(name)
                continue
            # make sure smarts string can be parsed
            try:
                self.smarts_parser.parse(smarts_string)
            except lark.ParseError as ex:
                if " col " in ex.args[0]:
                    column = ex.args[0][ex.args[0].find(" col ") + 5 :].strip()
                    column = " at character {} of {}".format(
                        column, smarts_string
                    )
                else:
                    column = ""
                malformed = "Malformed SMARTS string{} on line {}".format(column, entry.sourceline),
                errors.append(malformed)
                continue
            except lark.UnexpectedCharacters as e:
                now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
                m1="\t\tERROR: Unexpected Characters: \n"
                m2 ="\t\t{}".format(e)
                m3 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
                m = "\t\t" + len(m1) * "*" + "\n"
                print("\n" + m + m1 + m2 + m3 + m) if self._logger is None else \
                    self._logger.info("\n" + m + m1 + m2 + m3 + m)
                exit()

            # make sure referenced labels exist
            smarts_graph = SMARTSGraph(
                smarts_string,
                parser=self.smarts_parser,
                name=name,
                overrides=entry.attrib.get("overrides"),
            )
            for atom_expr in nx.get_node_attributes(smarts_graph, name="atom").values():
                labels = atom_expr.find_data("has_label")
                for label in labels:
                    atom_type = label.children[0][1:]
                    if atom_type not in self.atom_type_names:
                        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
                        m1 = "\t\tERROR: Reference to undefined atomtype '{}' in SMARTS string '{}'" \
                             " at line {}\n".format(atom_type, entry.attrib["def"], entry.sourceline)
                        m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
                        m = "\t\t" + len(m1) * "*" + "\n"
                        print("\n" + m + m1 + m2 + m) if self._logger is None \
                            else self._logger.info("\n" + m + m1 + m2 + m)
                        exit()

        if missing_smarts and debug:
            warn(
                "The following atom types do not have smarts definitions: {}".format(
                    ", ".join(missing_smarts)
                ),
                ValidationWarning,
            )
        if missing_smarts and not debug:

            now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
            m1 = "\t\tWARN: There are {} atom types that are missing a smarts definition. " \
                 "To view the missing atom types, re-run with debug=True when " \
                 "applying the forcefield.\n".format(len(missing_smarts))
            m = "\t\t" + len(m1) * "-" + "\n"
            print("\n" + m + m1 + m) if self._logger is None else self._logger.info("\n" + m + m1 + m)

    # =============================================================================================
    def validate_overrides(self):
        """Assert all overrides are defined elsewhere in force field."""
        errors = []
        for entry in self.atom_types:
            overrides = entry.attrib.get("overrides")
            if not overrides:
                continue
            overridden_types = [at.strip() for at in overrides.split(",") if at]
            for atom_type in overridden_types:
                if atom_type not in self.atom_type_names:
                    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
                    m1 = "\t\tERROR: Reference to undefined atomtype '{}' in 'overrides' '{}'" \
                         " at line {}\n".format(atom_type, entry.attrib["overrides"], entry.sourceline)
                    m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
                    m = "\t\t" + len(m1) * "*" + "\n"
                    print("\n" + m + m1 + m2 + m) if self._logger is None else self._logger.info("\n" + m + m1 + m2 + m)
                    exit()

    # =============================================================================================
    def get_all_forces_cj(self):

        """
        Get all labels in the XML file. This is useful to consider extra potential terms, i.e: PeriodicToxvaerdForce

        Args:

        Returns:
            A list of the tags for all Forces included in the XML

        """

        all_forces_included_xml = list()
        root = self.ff_tree.getroot()
        all_elements = self.ff_tree.getroot().findall("*")
        for item in all_elements:
            if item.tag.find("Force") != -1:
                all_forces_included_xml.append(item.tag)

        return all_forces_included_xml

