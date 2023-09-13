# =============================================================================
def print_header(version, logger_log=None):
    msg = """
    ***********************************************************************
              Replicate a molecule or polymer chain (ReMoPo)
                          top2forcefieldxml tool
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
        
        This tool converts a top/itp files produced with Q-Force toolkit
        (https://github.com/selimsami/qforce) to a xml forcefield to be used 
        with replicate_polymer tool. 

        This software is distributed under the terms of the
        GNU General Public License v3.0 (GNU GPLv3). A copy of
        the license (LICENSE.txt) is included with this distribution.

    ***********************************************************************
        """.format(version)

    print(msg) if logger_log is None else logger_log.info(msg)