import argparse
import argcomplete
import os
import datetime
import sys


# =============================================================================
def print_header(version, logger_log=None):
    msg = """
    ***********************************************************************
                 Remove hydrogens from the input file
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
    m += "\t\t\tremove_hydrogens".format(os.path.split(sys.argv[0])[1])
    m += m1+"\n"
    print(m) if logger is None else logger.info(m)


# =============================================================================
def parse_arguments():
    import time

    desc = """ Replicate and apply a force field to a polymer or molecule.\n"""

    parser = argparse.ArgumentParser(description=desc)
    # group1 = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument("-p", "--pdb", dest="pdbfile",
                        help="A pdb file containing the structure to be remove the hydrogens.",
                        action="store", metavar="INPUT_FILE", required=True)

    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    # Check for existing files:
    if not os.path.isfile(args.pdbfile):
        print(desc)
        time.sleep(.25)
        parser.error("PDB file must exist!!!!")

    return args
