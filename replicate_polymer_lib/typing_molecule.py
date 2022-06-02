import numpy as np
from replicate_polymer_lib.forcefield.Forcefield import Forcefield


def typing_molecule(filename, molecule, forcefield_xml_files, overwrite=True,
                    combining_rule="lorentz", type_kwargs=None, debug=False, logger=None, **kwargs):

    """

    Args:
        filename: Name of the file to write the topology using a Parmed instance
        molecule: A ReadPDBFormat object from topology library
        forcefield_xml_files: List of XML files containing force field parameters
        overwrite : bool, optional, default=False
                    Overwrite if the filename already exists
        combining_rule: Non-bonded combining rule
        type_kwargs:  dict, optional, default=None
                    Keyword arguments to provide to `Forcefield.apply`.
        debug: Debug FF info
        logger: Logger object

    Returns:

    """

    # Parse force field xml file
    ff = Forcefield(forcefield_files=forcefield_xml_files, debug=False, logger=logger)

    # Extra arguments
    if not type_kwargs:
        type_kwargs = {}

    # Non-periodic systems --> Create a box
    is_all_zero = np.all((molecule._unitcell == 0))
    if is_all_zero:
        xmin = min(molecule._universe.coord.positions[:, 0])
        ymin = min(molecule._universe.coord.positions[:, 1])
        zmin = min(molecule._universe.coord.positions[:, 2])
        xmax = max(molecule._universe.coord.positions[:, 0])
        ymax = max(molecule._universe.coord.positions[:, 1])
        zmax = max(molecule._universe.coord.positions[:, 2])
        molecule._unitcell = [[xmax-xmin + 5., 0.0, 0.0],
                              [0.0, ymax - ymin + 5., 0.0],
                              [0.0, 0.0, zmax - zmin + 5.]]
        molecule._boxlength = [xmax-xmin + 5., ymax-ymin + 5., zmax-zmin + 5.]
        molecule._boxangle = [90.0, 90.0, 90.0]


    pmd_structure = ff.apply(molecule, **type_kwargs)
    # pmd_structure.combining_rule = combining_rule

    total_charge = sum([atom.charge for atom in pmd_structure])
    if round(total_charge, 4) != 0.0:
        m1 = ("\t\tWARN: System is not charge neutral. Total charge is {}.\n".format(total_charge))
        m2 = "\t\tCheck charges!!!!!\n"
        ll = int(len(m1) / 1)
        m = "\t\t" + ll * "-" + "\n"
        print(m + m1 + m2 + m) if logger is None else logger.info(m + m1 + m2 + m)

    pmd_structure.save(filename, overwrite=overwrite, **kwargs)

