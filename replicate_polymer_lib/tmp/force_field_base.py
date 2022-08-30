import os
from collections import defaultdict
from copy import deepcopy
import xml.etree.ElementTree as ET
from replicate_polymer_lib.forcefield import unit


# =============================================================================================
def _convertparametertonumber(param):
    if unit.is_quantity(param):
        if param.unit.is_compatible(unit.bar):
            return param / unit.bar
        return param.value_in_unit_system(unit.md_unit_system)
    return float(param)


#==============================================================================================
# A map of functions to parse elements of the XML file.
parsers = {}


# ==============================================================================================
class FFBase():

    """A ForceField constructs a like-OpenMM System objects based on a Topology."""

    def __init__(self, *files):

        """Load one or more XML files and create a ForceField object based on them.

         Parameters
         ----------
         files : list
             A list of XML files defining the force field.  Each entry may
             be an absolute file path, a path relative to the current working
             directory, a path relative to this module's data subdirectory
             (for built in force fields), or an open file-like object with a
             read() method from which the forcefield XML data can be loaded.
         """
        self._atomTypes = {}
        self._templates = {}
        self._patches = {}
        self._templatePatches = {}
        self._templateSignatures = {None:[]}
        self._atomClasses = {'':set()}
        self._forces = []
        self._scripts = []
        self._templateMatchers = []
        self._templateGenerators = []
        self.load_file(files)

    # =============================================================================================
    def load_file(self, files, resname_prefix=''):

        """
        Load an XML file and add the definitions from it to this ForceField from OpenMM

        Args:
            files: string or file or tuple
                An XML file or tuple of XML files containing force field definitions.
                Each entry may be either an absolute file path, a path relative to the current working
                directory, a path relative to this module's data subdirectory (for
                built in force fields), or an open file-like object with a read()
                method from which the forcefield XML data can be loaded.
            resname_prefix : string
                An optional string to be prepended to each residue name found in the
                loaded files.
        """

        if isinstance(files, tuple):
            files = list(files)
        else:
            files = [files]

        trees = []

        i = 0
        while i < len(files):
            file = files[i]
            tree = None
            try:
                # this handles either filenames or open file-like objects
                tree = ET.parse(file)
            except Exception as e:
                # Fail with an error message about which file could not be read.
                # TODO: Also handle case where fallback to 'data' directory encounters problems,
                # but this is much less worrisome because we control those files.
                msg = str(e) + '\n'
                if hasattr(file, 'name'):
                    filename = file.name
                else:
                    filename = str(file)
                msg += "\t\tForceField.loadFile() encountered an error reading file '%s'\n" % filename
                print("\n\t\t"+msg) if self._logger is None else self._logger.info("\n\t\t"+msg)
                exit()
            if tree is None:
                m1 = '\t\tERROR: Could not locate file "{0:s}"\n'.format(file)
                m2 = "\t\tMolecule cannot be typed!!!!!. Exiting....\n"
                m = "\t\t" + len(m1) * "*" + "\n"
                print(m+m1+m2+m) if self._logger is None else self._logger.info(m+m1+m2+m)
                exit()
            trees.append(tree)
            i += 1

            # Process includes in this file.
            if isinstance(file, str):
                parent_dir = os.path.dirname(file)
            else:
                parent_dir = ''
            for included in tree.getroot().findall('Include'):
                include_file = included.attrib['file']
                joined = os.path.join(parent_dir, include_file)
                if os.path.isfile(joined):
                    include_file = joined
                if include_file not in files:
                    files.append(include_file)

        # Load the atom types.
        for tree in trees:
            if tree.getroot().find('AtomTypes') is not None:
                for itype in tree.getroot().find('AtomTypes').findall('Type'):
                    self.registerAtomType(itype.attrib)

        # Load the residue templates.
        for tree in trees:
            if tree.getroot().find('Residues') is not None:
                for residue in tree.getroot().find('Residues').findall('Residue'):
                    resName = resname_prefix+residue.attrib['name']
                    template = self._TemplateData(resName)
                    if 'override' in residue.attrib:
                        template.overrideLevel = int(residue.attrib['override'])
                    if 'rigidWater' in residue.attrib:
                        template.rigidWater = (residue.attrib['rigidWater'].lower() == 'true')
                    for key in residue.attrib:
                        template.attributes[key] = residue.attrib[key]
                    atomIndices = template.atomIndices
                    for ia, atom in enumerate(residue.findall('Atom')):
                        params = {}
                        for key in atom.attrib:
                            if key not in ('name', 'type'):
                                params[key] = _convertparametertonumber(atom.attrib[key])
                        atomName = atom.attrib['name']
                        if atomName in atomIndices:
                            raise ValueError('Residue '+resName+' contains multiple atoms named '+atomName)
                        typeName = atom.attrib['type']
                        atomIndices[atomName] = ia
                        template.atoms.append(self._TemplateAtomData(atomName, typeName, self._atomTypes[typeName].element, params))
                    for site in residue.findall('VirtualSite'):
                        template.virtualSites.append(self._VirtualSiteData(site, atomIndices))
                    for bond in residue.findall('Bond'):
                        if 'atomName1' in bond.attrib:
                            template.addBondByName(bond.attrib['atomName1'], bond.attrib['atomName2'])
                        else:
                            template.addBond(int(bond.attrib['from']), int(bond.attrib['to']))
                    for bond in residue.findall('ExternalBond'):
                        if 'atomName' in bond.attrib:
                            template.addExternalBondByName(bond.attrib['atomName'])
                        else:
                            template.addExternalBond(int(bond.attrib['from']))
                    for patch in residue.findall('AllowPatch'):
                        patchName = patch.attrib['name']
                        if ':' in patchName:
                            colonIndex = patchName.find(':')
                            self.registerTemplatePatch(resName, patchName[:colonIndex], int(patchName[colonIndex+1:])-1)
                        else:
                            self.registerTemplatePatch(resName, patchName, 0)
                    self.registerResidueTemplate(template)

        # Load the patch definitions.
        for tree in trees:
            if tree.getroot().find('Patches') is not None:
                for patch in tree.getroot().find('Patches').findall('Patch'):
                    patchName = patch.attrib['name']
                    if 'residues' in patch.attrib:
                        numResidues = int(patch.attrib['residues'])
                    else:
                        numResidues = 1
                    patchData = self._PatchData(patchName, numResidues)
                    for key in patch.attrib:
                        patchData.attributes[key] = patch.attrib[key]
                    for atom in patch.findall('AddAtom'):
                        params = {}
                        for key in atom.attrib:
                            if key not in ('name', 'type'):
                                params[key] = _convertparametertonumber(atom.attrib[key])
                        atomName = atom.attrib['name']
                        if atomName in patchData.allAtomNames:
                            raise ValueError('Patch '+patchName+' contains multiple atoms named '+atomName)
                        patchData.allAtomNames.add(atomName)
                        atomDescription = self._PatchAtomData(atomName)
                        typeName = atom.attrib['type']
                        patchData.addedAtoms[atomDescription.residue].append(self._TemplateAtomData(atomDescription.name, typeName, self._atomTypes[typeName].element, params))
                    for atom in patch.findall('ChangeAtom'):
                        params = {}
                        for key in atom.attrib:
                            if key not in ('name', 'type'):
                                params[key] = _convertparametertonumber(atom.attrib[key])
                        atomName = atom.attrib['name']
                        if atomName in patchData.allAtomNames:
                            raise ValueError('Patch '+patchName+' contains multiple atoms named '+atomName)
                        patchData.allAtomNames.add(atomName)
                        atomDescription = self._PatchAtomData(atomName)
                        typeName = atom.attrib['type']
                        patchData.changedAtoms[atomDescription.residue].append(self._TemplateAtomData(atomDescription.name, typeName, self._atomTypes[typeName].element, params))
                    for atom in patch.findall('RemoveAtom'):
                        atomName = atom.attrib['name']
                        if atomName in patchData.allAtomNames:
                            raise ValueError('Patch '+patchName+' contains multiple atoms named '+atomName)
                        patchData.allAtomNames.add(atomName)
                        atomDescription = self._PatchAtomData(atomName)
                        patchData.deletedAtoms.append(atomDescription)
                    for bond in patch.findall('AddBond'):
                        atom1 = self._PatchAtomData(bond.attrib['atomName1'])
                        atom2 = self._PatchAtomData(bond.attrib['atomName2'])
                        patchData.addedBonds.append((atom1, atom2))
                    for bond in patch.findall('RemoveBond'):
                        atom1 = self._PatchAtomData(bond.attrib['atomName1'])
                        atom2 = self._PatchAtomData(bond.attrib['atomName2'])
                        patchData.deletedBonds.append((atom1, atom2))
                    for bond in patch.findall('AddExternalBond'):
                        atom = self._PatchAtomData(bond.attrib['atomName'])
                        patchData.addedExternalBonds.append(atom)
                    for bond in patch.findall('RemoveExternalBond'):
                        atom = self._PatchAtomData(bond.attrib['atomName'])
                        patchData.deletedExternalBonds.append(atom)
                    atomIndices = dict((atom.name, i) for i, atom in enumerate(patchData.addedAtoms[atomDescription.residue]+patchData.changedAtoms[atomDescription.residue]))
                    for site in patch.findall('VirtualSite'):
                        patchData.virtualSites[atomDescription.residue].append(self._VirtualSiteData(site, atomIndices))
                    for residue in patch.findall('ApplyToResidue'):
                        name = residue.attrib['name']
                        if ':' in name:
                            colonIndex = name.find(':')
                            self.registerTemplatePatch(name[colonIndex+1:], patchName, int(name[:colonIndex])-1)
                        else:
                            self.registerTemplatePatch(name, patchName, 0)
                    self.registerPatch(patchData)

        # Load force definitions
        for tree in trees:
            for child in tree.getroot():
                if child.tag in parsers:
                    parsers[child.tag](child, self)

        # Load scripts
        for tree in trees:
            for node in tree.getroot().findall('Script'):
                self.registerScript(node.text)

        # Execute initialization scripts.
        for tree in trees:
            for node in tree.getroot().findall('InitializationScript'):
                exec(node.text, locals())

    # =============================================================================================
    def registerGenerator(self, generator):
        """Register a new generator."""
        self._forces.append(generator)

    # =============================================================================================
    def registerAtomType(self, parameters):
        """Register a new atom type."""
        name = parameters["name"]
        if name in self._atomTypes:
            raise ValueError(
                "Found multiple definitions for atom type: " + name
            )
        atom_class = parameters["class"]
        mass = _convertparametertonumber(parameters["mass"])
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
    def registerTemplatePatch(self, residue, patch, patchResidueIndex):
        """Register that a particular patch can be used with a particular residue."""
        if residue not in self._templatePatches:
            self._templatePatches[residue] = set()
        self._templatePatches[residue].add((patch, patchResidueIndex))

    # =============================================================================================
    def registerResidueTemplate(self, template):
        """Register a new residue template."""
        if template.name in self._templates:
            # There is already a template with this name, so check the override levels.

            existingTemplate = self._templates[template.name]
            if template.overrideLevel < existingTemplate.overrideLevel:
                # The existing one takes precedence, so just return.
                return
            if template.overrideLevel > existingTemplate.overrideLevel:
                # We need to delete the existing template.
                del self._templates[template.name]
                existingSignature = _create_residue_signature([atom.element for atom in existingTemplate.atoms])
                self._templateSignatures[existingSignature].remove(existingTemplate)
            else:
                raise ValueError('Residue template %s with the same override level %d already exists.' % (template.name, template.overrideLevel))

        # Register the template.

        self._templates[template.name] = template
        signature = _create_residue_signature([atom.element for atom in template.atoms])
        if signature in self._templateSignatures:
            self._templateSignatures[signature].append(template)
        else:
            self._templateSignatures[signature] = [template]

    # =============================================================================================
    def _findAtomTypes(self, attrib, num):
        """Parse the attributes on an XML tag to find the set of atom types for each atom it involves.

        Parameters
        ----------
        attrib : dict of attributes
            The dictionary of attributes for an XML parameter tag.
        num : int
            The number of atom specifiers (e.g. 'class1' through 'class4') to extract.

        Returns
        -------
        types : list
            A list of atom types that match.

        """
        types = []
        for i in range(num):
            if num == 1:
                suffix = ''
            else:
                suffix = str(i+1)
            classAttrib = 'class'+suffix
            typeAttrib = 'type'+suffix
            if classAttrib in attrib:
                if typeAttrib in attrib:
                    raise ValueError('Specified both a type and a class for the same atom: '+str(attrib))
                if attrib[classAttrib] not in self._atomClasses:
                    types.append(None) # Unknown atom class
                else:
                    types.append(self._atomClasses[attrib[classAttrib]])
            elif typeAttrib in attrib:
                if attrib[typeAttrib] == '':
                    types.append(self._atomClasses[''])
                elif attrib[typeAttrib] not in self._atomTypes:
                    types.append(None) # Unknown atom type
                else:
                    types.append([attrib[typeAttrib]])
            else:
                types.append(None) # Unknown atom type
        return types

    # =============================================================================================
    class _TemplateData(object):
        """Inner class used to encapsulate data about a residue template definition."""
        def __init__(self, name):
            self.name = name
            self.atoms = []
            self.atomIndices = {}
            self.virtualSites = []
            self.bonds = []
            self.externalBonds = []
            self.overrideLevel = 0
            self.rigidWater = True
            self.attributes = {}

        def getAtomIndexByName(self, atom_name):
            """Look up an atom index by atom name, providing a helpful error message if not found."""
            index = self.atomIndices.get(atom_name, None)
            if index is not None:
                return index

            # Provide a helpful error message if atom name not found.
            msg =  "Atom name '%s' not found in residue template '%s'.\n" % (atom_name, self.name)
            msg += "Possible atom names are: %s" % str(list(map(lambda x: x.name, self.atoms)))
            raise ValueError(msg)

        def addAtom(self, atom):
            self.atoms.append(atom)
            self.atomIndices[atom.name] = len(self.atoms)-1

        def addBond(self, atom1, atom2):
            """Add a bond between two atoms in a template given their indices in the template."""
            self.bonds.append((atom1, atom2))
            self.atoms[atom1].bondedTo.append(atom2)
            self.atoms[atom2].bondedTo.append(atom1)

        def addBondByName(self, atom1_name, atom2_name):
            """Add a bond between two atoms in a template given their atom names."""
            atom1 = self.getAtomIndexByName(atom1_name)
            atom2 = self.getAtomIndexByName(atom2_name)
            self.addBond(atom1, atom2)

        def addExternalBond(self, atom_index):
            """Designate that an atom in a residue template has an external bond, using atom index within template."""
            self.externalBonds.append(atom_index)
            self.atoms[atom_index].externalBonds += 1

        def addExternalBondByName(self, atom_name):
            """Designate that an atom in a residue template has an external bond, using atom name within template."""
            atom = self.getAtomIndexByName(atom_name)
            self.addExternalBond(atom)

        def areParametersIdentical(self, template2, matchingAtoms, matchingAtoms2):
            """Get whether this template and another one both assign identical atom types and parameters to all atoms.

            Parameters
            ----------
            template2: _TemplateData
                the template to compare this one to
            matchingAtoms: list
                the indices of atoms in this template that atoms of the residue are matched to
            matchingAtoms2: list
                the indices of atoms in template2 that atoms of the residue are matched to
            """
            atoms1 = [self.atoms[m] for m in matchingAtoms]
            atoms2 = [template2.atoms[m] for m in matchingAtoms2]
            if any(a1.type != a2.type or a1.parameters != a2.parameters for a1,a2 in zip(atoms1, atoms2)):
                return False
            # Properly comparing virtual sites really needs a much more complicated analysis.  This simple check
            # could easily fail for templates containing vsites, even if they're actually identical.  Since we
            # currently have no force fields that include both patches and vsites, I'm not going to worry about it now.
            if self.virtualSites != template2.virtualSites:
                return False
            return True

    # =============================================================================================
    class _TemplateAtomData(object):
        """Inner class used to encapsulate data about an atom in a residue template definition."""
        def __init__(self, name, itype, element, parameters={}):
            self.name = name
            self.type = itype
            self.element = element
            self.parameters = parameters
            self.bondedTo = []
            self.externalBonds = 0

    # =============================================================================================
    class _VirtualSiteData(object):
        """Inner class used to encapsulate data about a virtual site."""
        def __init__(self, node, atomIndices):
            attrib = node.attrib
            self.type = attrib['type']
            if self.type == 'average2':
                numAtoms = 2
                self.weights = [float(attrib['weight1']), float(attrib['weight2'])]
            elif self.type == 'average3':
                numAtoms = 3
                self.weights = [float(attrib['weight1']), float(attrib['weight2']), float(attrib['weight3'])]
            elif self.type == 'outOfPlane':
                numAtoms = 3
                self.weights = [float(attrib['weight12']), float(attrib['weight13']), float(attrib['weightCross'])]
            elif self.type == 'localCoords':
                numAtoms = 0
                self.originWeights = []
                self.xWeights = []
                self.yWeights = []
                while ('wo%d' % (numAtoms+1)) in attrib:
                    numAtoms += 1
                    self.originWeights.append(float(attrib['wo%d' % numAtoms]))
                    self.xWeights.append(float(attrib['wx%d' % numAtoms]))
                    self.yWeights.append(float(attrib['wy%d' % numAtoms]))
                self.localPos = [float(attrib['p1']), float(attrib['p2']), float(attrib['p3'])]
            else:
                raise ValueError('Unknown virtual site type: %s' % self.type)
            if 'siteName' in attrib:
                self.index = atomIndices[attrib['siteName']]
                self.atoms = [atomIndices[attrib['atomName%d'%(i+1)]] for i in range(numAtoms)]
            else:
                self.index = int(attrib['index'])
                self.atoms = [int(attrib['atom%d'%(i+1)]) for i in range(numAtoms)]
            if 'excludeWith' in attrib:
                self.excludeWith = int(attrib['excludeWith'])
            else:
                self.excludeWith = self.atoms[0]

    # =============================================================================================
    class _PatchData(object):
        """Inner class used to encapsulate data about a patch definition."""
        def __init__(self, name, numResidues):
            self.name = name
            self.numResidues = numResidues
            self.addedAtoms = [[] for i in range(numResidues)]
            self.changedAtoms = [[] for i in range(numResidues)]
            self.deletedAtoms = []
            self.addedBonds = []
            self.deletedBonds = []
            self.addedExternalBonds = []
            self.deletedExternalBonds = []
            self.allAtomNames = set()
            self.virtualSites = [[] for i in range(numResidues)]
            self.attributes = {}

        def createPatchedTemplates(self, templates):
            """Apply this patch to a set of templates, creating new modified ones."""
            if len(templates) != self.numResidues:
                raise ValueError("Patch '%s' expected %d templates, received %d", (self.name, self.numResidues, len(templates)))

            # Construct a new version of each template.
            newTemplates = []
            for index, template in enumerate(templates):
                newTemplate = self._TemplateData("%s-%s" % (template.name, self.name))
                newTemplates.append(newTemplate)

                # Build the list of atoms in it.

                for atom in template.atoms:
                    if not any(deleted.name == atom.name and deleted.residue == index for deleted in self.deletedAtoms):
                        newTemplate.addAtom(self._TemplateAtomData(atom.name, atom.type, atom.element, atom.parameters))
                for atom in self.addedAtoms[index]:
                    if any(a.name == atom.name for a in newTemplate.atoms):
                        raise ValueError("Patch '%s' adds an atom with the same name as an existing atom: %s" % (self.name, atom.name))
                    newTemplate.addAtom(self._TemplateAtomData(atom.name, atom.type, atom.element, atom.parameters))
                oldAtomIndex = dict([(atom.name, i) for i, atom in enumerate(template.atoms)])
                newAtomIndex = dict([(atom.name, i) for i, atom in enumerate(newTemplate.atoms)])
                for atom in self.changedAtoms[index]:
                    if atom.name not in newAtomIndex:
                        raise ValueError("Patch '%s' modifies nonexistent atom '%s' in template '%s'" % (self.name, atom.name, template.name))
                    newTemplate.atoms[newAtomIndex[atom.name]] = self._TemplateAtomData(atom.name, atom.type, atom.element, atom.parameters)

                # Copy over the virtual sites, translating the atom indices.
                indexMap = dict([(oldAtomIndex[name], newAtomIndex[name]) for name in newAtomIndex if name in oldAtomIndex])
                for site in template.virtualSites:
                    if site.index in indexMap and all(i in indexMap for i in site.atoms):
                        newSite = deepcopy(site)
                        newSite.index = indexMap[site.index]
                        newSite.atoms = [indexMap[i] for i in site.atoms]
                        newTemplate.virtualSites.append(newSite)

                # Build the lists of bonds and external bonds.

                atomMap = dict([(template.atoms[i], indexMap[i]) for i in indexMap])
                deletedBonds = [(atom1.name, atom2.name) for atom1, atom2 in self.deletedBonds if atom1.residue == index and atom2.residue == index]
                for atom1, atom2 in template.bonds:
                    a1 = template.atoms[atom1]
                    a2 = template.atoms[atom2]
                    if a1 in atomMap and a2 in atomMap and (a1.name, a2.name) not in deletedBonds and (a2.name, a1.name) not in deletedBonds:
                        newTemplate.addBond(atomMap[a1], atomMap[a2])
                deletedExternalBonds = [atom.name for atom in self.deletedExternalBonds if atom.residue == index]
                for atom in template.externalBonds:
                    if template.atoms[atom].name not in deletedExternalBonds:
                        newTemplate.addExternalBond(indexMap[atom])
                for atom1, atom2 in self.addedBonds:
                    if atom1.residue == index and atom2.residue == index:
                        newTemplate.addBondByName(atom1.name, atom2.name)
                    elif atom1.residue == index:
                        newTemplate.addExternalBondByName(atom1.name)
                    elif atom2.residue == index:
                        newTemplate.addExternalBondByName(atom2.name)
                for atom in self.addedExternalBonds:
                    newTemplate.addExternalBondByName(atom.name)

                # Add new virtual sites.

                indexMap = dict((i, newAtomIndex[atom.name]) for i, atom in enumerate(self.addedAtoms[index]+self.changedAtoms[index]))
                for site in self.virtualSites[index]:
                    newSite = deepcopy(site)
                    newSite.index = indexMap[site.index]
                    newSite.atoms = [indexMap[i] for i in site.atoms]
                    newSite.excludeWith = indexMap[site.excludeWith]
                    newTemplate.virtualSites = [site for site in newTemplate.virtualSites if site.index != newSite.index]
                    newTemplate.virtualSites.append(newSite)
            return newTemplates

    # =============================================================================================
    class _PatchAtomData(object):
        """Inner class used to encapsulate data about an atom in a patch definition."""
        def __init__(self, description):
            if ':' in description:
                colonIndex = description.find(':')
                self.residue = int(description[:colonIndex])-1
                self.name = description[colonIndex+1:]
            else:
                self.residue = 0
                self.name = description


# =============================================================================================
def _create_residue_signature(elements):
    """Create a signature for a residue based on the elements of the atoms it contains."""
    counts = _count_residue_atoms(elements)
    sig = []
    for c in counts:
        if c is not None:
            sig.append((c, counts[c]))
    sig.sort(key=lambda x: -x[0].mass)

    # Convert it to a string.

    s = ''
    for element, count in sig:
        s += element.symbol+str(count)
    return s


# =============================================================================================
def _count_residue_atoms(elements):
    """Count the number of atoms of each element in a residue."""
    counts = {}
    for element in elements:
        if element in counts:
            counts[element] += 1
        else:
            counts[element] = 1
    return counts

# The following classes are generators that know how to create Force subclasses and add them to a System that is being
# created.  Each generator class must define two methods: 1) a static method that takes an etree Element and a ForceField,
# and returns the corresponding generator object; 2) a createForce() method that constructs the Force object and adds it
# to the System.  The static method should be added to the parsers map.

class HarmonicBondGenerator(object):
    """A HarmonicBondGenerator constructs a HarmonicBondForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.bondsForAtomType = defaultdict(set)
        self.types1 = []
        self.types2 = []
        self.length = []
        self.k = []

    def registerBond(self, parameters):
        types = self.ff._findAtomTypes(parameters, 2)
        if None not in types:
            index = len(self.types1)
            self.types1.append(types[0])
            self.types2.append(types[1])
            for t in types[0]:
                self.bondsForAtomType[t].add(index)
            for t in types[1]:
                self.bondsForAtomType[t].add(index)
            self.length.append(_convertparametertonumber(parameters['length']))
            self.k.append(_convertparametertonumber(parameters['k']))

    @staticmethod
    def parseElement(element, ff):
        existing = [f for f in ff._forces if isinstance(f, HarmonicBondGenerator)]
        if len(existing) == 0:
            generator = HarmonicBondGenerator(ff)
            ff.registerGenerator(generator)
        else:
            generator = existing[0]
        for bond in element.findall('Bond'):
            generator.registerBond(bond.attrib)

    # def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
    #     existing = [f for f in sys.getForces() if type(f) == mm.HarmonicBondForce]
    #     if len(existing) == 0:
    #         force = mm.HarmonicBondForce()
    #         sys.addForce(force)
    #     else:
    #         force = existing[0]
    #     for bond in data.bonds:
    #         type1 = data.atomType[data.atoms[bond.atom1]]
    #         type2 = data.atomType[data.atoms[bond.atom2]]
    #         for i in self.bondsForAtomType[type1]:
    #             types1 = self.types1[i]
    #             types2 = self.types2[i]
    #             if (type1 in types1 and type2 in types2) or (type1 in types2 and type2 in types1):
    #                 bond.length = self.length[i]
    #                 if bond.isConstrained:
    #                     data.addConstraint(sys, bond.atom1, bond.atom2, self.length[i])
    #                 if self.k[i] != 0:
    #                     # flexibleConstraints allows us to add parameters even if the DOF is
    #                     # constrained
    #                     if not bond.isConstrained or args.get('flexibleConstraints', False):
    #                         force.addBond(bond.atom1, bond.atom2, self.length[i], self.k[i])
    #                 break

parsers["HarmonicBondForce"] = HarmonicBondGenerator.parseElement
