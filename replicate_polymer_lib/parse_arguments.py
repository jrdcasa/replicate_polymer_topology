import argparse
import argcomplete
import os
import sys


# =============================================================================
def print_header(version, logger_log=None):
    msg = """
    ***********************************************************************
              Replicate a molecule or polymer chain (ReMoPo)
              ----------------------------------------------

                                Version {}

                              Dr. Javier Ramos
                      Macromolecular Physics Department
                Instituto de Estructura de la Materia (IEM-CSIC)
                               Madrid (Spain)

        ReMoPo is an open-source python library to quickly replicate a molecule
        from a pdb file containing the seed molecule. After that, the script assigns
        the force field parameters using the foyer library (https://mosdef.org/)
        This program generates GROMACS and LAMMPS files to run simulations.

        This software is distributed under the terms of the
        GNU General Public License v3.0 (GNU GPLv3). A copy of
        the license (LICENSE.txt) is included with this distribution.

    ***********************************************************************
        """.format(version)

    print(msg) if logger_log is None else logger_log.info(msg)


# =============================================================================
def command_info(opts, logger=None):

    m1 = ""
    for item in sys.argv[1:]:
        m1 += " {}".format(item)
    m = "\n\t\tCommand line: \n"
    m += "\t\t\tpython {}".format(os.path.split(sys.argv[0])[1])
    m += m1+"\n"
    m += "\t\t\t         or\n"
    m += "\t\t\treplicate_polymer".format(os.path.split(sys.argv[0])[1])
    m += m1+"\n"
    print(m) if logger is None else logger.info(m)

    if opts.verbose:
        m1 = "\t\tWARN: Verbose checking of angles and dihedral is activated. " \
             "\n\t\tThis is time consuming. Perform only if there is angles or dihedrals without assignment.\n"
        ll = int(len(m1)/1.5)
        m = "\t\t" + ll * "-" + "\n"
        print(m+m1+m) if logger is None else logger.info(m+m1+m)


# =============================================================================
def parse_arguments():
    import time

    desc = """ Replicate and apply a force field to a polymer or molecule.\n"""

    parser = argparse.ArgumentParser(description=desc)
    # group1 = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument("-p", "--pdb", dest="pdbfile",
                        help="A pdb file containing the structure to be replicated.",
                        action="store", metavar="PDB_FILE", required=True)
    # group2 = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument("-f", "--forcefield", dest="xmlfile",
                        help="A XML file containing the force field parameters using the FOYER format",
                        action="store", metavar="XML_FILE", required=True)

    parser.add_argument("--images", dest="images", nargs=3,
                        help="A list with the number of images to replicate in the x, y and z axis.",
                        action="store", metavar=('image_x', 'image_y', 'image_z'), required=True, default=None)

    parser.add_argument("-e", "--engine", dest="engine",
                        help="MD package to perform the calculations.",
                        action="store", metavar="MDENGINE", required=False, default="gromacs")

    # group2 = parser.add_mutually_exclusive_group(required=False)
    parser.add_argument("--noh", dest="hydrogens",
                        help="Remove hydrogens for a united atom representation.",
                        action="store_true", required=False)

    parser.add_argument("--index", dest="index",
                        help="Indices of atoms to be removed from the PDB.",
                        action="store", metavar="INDEX", required=False, default=None)

    parser.add_argument("--boxlength", dest="boxlength", nargs=3,
                        help="Box lengths in nanometers.", type=float,
                        action="store", metavar=('a', 'b', 'c'), required=False, default=None)

    parser.add_argument("--boxangle", dest="boxangle", nargs=3,
                        help="Box angles in degrees.", type=float,
                        action="store", metavar=('alpha', 'beta', 'gamma'), required=False, default=None)

    parser.add_argument("--impropers", dest="impropers",
                        help="Write impropers to top file. The argument can be"
                             "\'Guess\' or the name of a file with the improper angles defined. "
                             "If None improper angles neither be guessed or read from a file",
                        action="store", required=False, default=None)

    parser.add_argument("--npairs", dest="npairs",
                        help="Monomer or residue inclusions. A value of 1 indicates that only non-bonded interactions"
                             "between the current residue (i) and the nearest-neighbours are "
                             "taken into account (i-1, i and i+1)", type=int,
                        action="store", required=False, default=None)

    parser.add_argument("--verbose", dest="verbose",
                        help="Verbose checking of angles and dihedral. This is time consuming. "
                             "Perform only if there is angles without assignment.",
                        action="store_true", required=False)

    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    # Check for existing files:
    if not os.path.isfile(args.pdbfile):
        print(desc)
        time.sleep(.25)
        parser.error("PDB file must exist!!!!")
    if not os.path.isfile(args.xmlfile):
        print(desc)
        time.sleep(.25)
        parser.error("XML file must exist!!!!")
    if not args.images:
        args.images = None
    else:
        # Check integers
        try:
            args.images = [int(x) for x in args.images]
        except ValueError:
            print(desc)
            print("Actual images list: ", args.images)
            time.sleep(.25)
            parser.error("Images must be a list of integers!!!!")

    if not args.boxlength:
        args.boxlength = None
    if not args.boxangle:
        args.boxangle = None
    allow_engines = ["gromacs", "lammps", "lammps_mc"]

    if args.engine.lower() not in allow_engines:
        print(desc)
        time.sleep(.25)
        line = "Allowed MD engines: " + ", ".join(allow_engines)
        parser.error(line)
    else:
        args.engine = args.engine.lower()

    if args.impropers is None:
        args.impropers = "No"
    elif args.impropers.upper() == "GUESS":
        args.impropers = "GUESS"
    else:
        if not os.path.isfile(args.impropers):
            line = "\nImproper file must exist. The format should follow the ndx format of GROMACS\n" \
                   "The file must start with the label [ impropers ]\n" \
                   "Each line is a improper angle.\n" \
                   "File {} is not found".format(args.impropers)
            parser.error(line)

    return args
