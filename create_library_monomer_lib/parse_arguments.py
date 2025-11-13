import argparse
import argcomplete
import os
import sys


# =============================================================================
def print_header(version, logger_log=None):
    msg = """
    ***********************************************************************
              Create itp monomers and apply a force field 
                   to be included in a monomer library
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
    m += '\t\t\t{}'.format(os.path.split(sys.argv[0])[1])
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

    desc = """ Create itp monomers and apply a force field to be included in a monomer library.\n"""

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("-p", "--pdb", dest="pdbfile", nargs="+",
                        help="A pdb file containing the structure to be replicated.",
                        action="store", metavar="PDB_FILE", required=True, default=None)

    parser.add_argument("-f", "--forcefield", dest="forcefield",
                        help="A directory or a file that stores the initial XML, ITP, or TOP files "
                             "containing the force field parameters in FOYER and GROMACS formats.",
                        action="store", metavar="FF_FILE", required=True)

    parser.add_argument("--pattern", dest="pattern",
                        help="A naming pattern for the directory to store the ff files for polyply.",
                        action="store", metavar="PATTERN", required=True)

    # group2 = parser.add_mutually_exclusive_group(required=False)
    parser.add_argument("--noh", dest="hydrogens",
                        help="Remove hydrogens for a united atom representation.",
                        action="store_true", required=False)

    parser.add_argument("--nodeleteff", dest="nodeleteff",
                        help="Delete all previoulsy *.ff files in the working directory.",
                        action="store_true", required=False, default=True)

    parser.add_argument("--debug", dest="debug",
                        help="Debug mode.",
                        action="store_true", required=False, default=False)

    parser.add_argument("--verbose", dest="verbose",
                        help="Verbose checking of angles and dihedral. This is time consuming. "
                             "Perform only if there is angles without assignment.",
                        action="store_true", required=False)

    parser.add_argument("-s", "--seq", dest="sequences", nargs="+",
                        help="Sequences to build the polymer chain. "
                             "The format follow is that defined in polyply",
                        action="store", metavar="SEQ_STR", required=False, default=None)

    parser.add_argument("-n", "--names", dest="namesystems", nargs="+",
                        help="Name of the systems to be created. "
                             "In case of you do not provide it, the name will be built using the provided PDB files.",
                        action="store", metavar="NAMES_STR", required=False, default=None)

    parser.add_argument("--nmols", dest="nmols", nargs="+",
                        help="Number of molecules for each system. ",
                        action="store", metavar="NMOLS", required=False, default=None)

    parser.add_argument("--gen_coords_opts", dest="gen_coords_opts", nargs="+",
                        help="Options for polyply gen_coords <density kg/cm3>, "
                             "<mir: max number of trys to place a residue>, "
                             "<mf: max force to allow in residue placing (kJ/mol*nm)>,"
                             "<nr: Number of residues to trace back when RW fails in first attempt> ",
                        action="store", metavar="PARAMS", required=False, default=[650, 1000, 40000, 5])

    parser.add_argument("--gmxpath", dest="gmxpath",
                        help="Path to the gmx program (GROMACS).",
                        action="store", metavar="GMX PATH", required=False, default=None)

    parser.add_argument("--charges", dest="charges",
                        help="A file containing the atomistic charges of each atom. "
                             "The order of the atoms must be the same that teh order in the pdb.",
                        action="store", metavar="GMX PATH", required=False, default=None)

    parser.add_argument("--residues", dest="residues",
                        help="A file containing the residue to which each atom belongs. "
                             "The order of the atoms must be the same as the order in the PDB",
                        action="store", metavar="GMX PATH", required=False, default=None)

    parser.add_argument("--noavgcharges", dest="noavgcharges",
                        help="Do not perform charge averaging on residues of the same type.",
                        action="store_true", required=False)

    parser.add_argument("--zerochargesres", dest="zerochargesres",
                        help="Remove excess of charge in the residues.",
                        action="store_true", required=False)

    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    # Check for existing files:
    for ipdb in args.pdbfile:
        if not os.path.isfile(ipdb):
            parser.error("\nERROR: File {} does not exist\n".format(ipdb))

    if args.sequences is not None:
        if len(args.sequences) != len(args.pdbfile):
            parser.error("\nERROR: Sequences length (--seq) must be equal to pdb length (--pdb)\n"
                         "  Seq length = {0:d}\n  PDB length = {1:d}".format(len(args.sequences), len(args.pdbfile)))

    if args.namesystems is not None:
        if len(args.namesystems) != len(args.pdbfile):
            parser.error("\nERROR: Name systems length (--namesystems) must be equal to pdb length (--pdb)\n"
                         "  Name system length = {0:d}\n  PDB length = {1:d}".
                         format(len(args.namesystems), len(args.pdbfile)))

    if args.nmols is not None:
        if len(args.nmols) != len(args.pdbfile):
            parser.error("\nERROR: Number of molecules length (--namesystems) must be equal to pdb length (--pdb)\n"
                         "  Name system length = {0:d}\n  PDB length = {1:d}".
                         format(len(args.nmols), len(args.pdbfile)))

    if not os.path.isdir(args.forcefield) and not os.path.isfile(args.forcefield):
        print(desc)
        time.sleep(.25)
        parser.error("Forcefield directory with XML, ITP or TOP files must exist!!!!")

    return args
