import argparse
import argcomplete
import os
import datetime
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
    m += "\t\t\tcreate_lammps_from_gromacs".format(os.path.split(sys.argv[0])[1])
    m += m1+"\n"
    print(m) if logger is None else logger.info(m)


# =============================================================================
def parse_arguments():
    import time

    desc = """ Create templates files for LAMMPS from GROMACS top and gro files.\n"""

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-g", "--gro", dest="grofile",
                        help="A gro file containing the coordinates of the system.",
                        action="store", metavar="GRO_FILE", required=True)
    parser.add_argument("-t", "--top", dest="topfile",
                        help="A top file containing the coordinates of the system.",
                        action="store", metavar="XML_FILE", required=True)
    parser.add_argument("--mc", dest="mc",
                        help="Lammps format for MC single chain program.",
                        action="store_true", required=False)
    parser.add_argument("-f", "--forcefield", dest="xmlfile",
                        help="A XML file containing the force field parameters using the FOYER format",
                        action="store", metavar="XML_FILE", required=True)

    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    # Check for existing files:
    if not os.path.isfile(args.grofile):
        print(desc)
        time.sleep(.25)
        parser.error("GRO file must exist!!!!")
    if not os.path.isfile(args.topfile):
        print(desc)
        time.sleep(.25)
        parser.error("TOP file must exist!!!!")
    if not os.path.isfile(args.xmlfile):
        print(desc)
        time.sleep(.25)
        parser.error("XML file must exist!!!!")

    return args
