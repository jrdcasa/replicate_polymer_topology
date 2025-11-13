import datetime
from remove_hydrogens_lib.parse_arguments import parse_arguments, print_header, command_info
from replicate_polymer_lib.logger import init_logger
from replicate_polymer_lib.check_connect_pdb import check_conect_pdb, check_resname_pdb, check_and_remove_ter_labels
from replicate_polymer_lib.remove_hydrogens import remove_hydrogens


# =============================================================================
def main_app(version):

    # Init logger
    logger = init_logger("Output", fileoutput="Info.log", append=False, inscreen=True)
    print_header(version, logger)

    # Parse arguments in the command line
    opts = parse_arguments()

    # Info
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    starttime = datetime.datetime.now()
    m = "\t\tStart Job at {} ============".format(now)
    print(m) if logger is None else logger.info(m)

    # Print command info
    command_info(opts, logger=logger)

    # Read the original pdb and check PDB Connect section in the PDB file
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

    m = "\n\t\tCheck resname and connect in the input file. ({})... ({})".format(opts.pdbfile, now)
    m += "\n\t\t{}".format(len(m)*"*")
    print(m) if logger is None else logger.info(m)
    filenamepdb = check_resname_pdb(opts.pdbfile)
    filenamepdb, bond_list_more100K = check_conect_pdb(filenamepdb, logger=logger)
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\tEnd Check resname and connect in the input file. ({})... ({})".format(opts.pdbfile, now)
    m += "\n\t\t{}".format(len(m)*"*")
    print(m) if logger is None else logger.info(m)
    #
    # Remove hydrogen atoms
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\tRemoving hydrogens in the pdb file of the seed molecules.({})... ({})".format(filenamepdb, now)
    print(m) if logger is None else logger.info(m)
    filenamepdb_new, filenamegro_new = remove_hydrogens(filenamepdb, bond_list_more100K)
    filenamepdb = filenamepdb_new

    filenamepdb = check_and_remove_ter_labels(filenamepdb, logger=logger)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    endtime = datetime.datetime.now()
    m = "\t\tJob  Done at {} ============\n".format(now)
    print(m) if logger is None else logger.info(m)
    m = "\t\tTotal time: {0:.2f} seconds".format((endtime-starttime).total_seconds())
    print(m) if logger is None else logger.info(m)


# =============================================================================
if __name__ == "__main__":

    __version__ = "1.1"
    main_app(__version__)
