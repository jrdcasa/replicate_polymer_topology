import datetime
from replicate_polymer_lib.logger import init_logger
from replicate_polymer_lib.prepare_lammps import clean_lammps
from create_lammps_from_gromacs_lib.parse_arguments import parse_arguments, print_header, command_info


# =============================================================================
def main_app(version):

    # Init logger
    logger = init_logger("Output", fileoutput="Info_lammps.log", append=False, inscreen=True)
    print_header(version, logger)
    #
    # Parse arguments in the command line
    opts = parse_arguments()
    #
    # Info
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    starttime = datetime.datetime.now()
    m = "\t\tStart Job at {} ============".format(now)
    print(m) if logger is None else logger.info(m)

    # Print command info
    command_info(opts, logger=logger)

    basefile_gro = opts.grofile
    basefile_top = opts.topfile
    if opts.mc:
        clean_lammps(basefile_gro, basefile_top, opts.xmlfile, create_lmp_mc_singlechain=True)
    else:
        clean_lammps(basefile_gro, basefile_top, opts.xmlfile, create_lmp_mc_singlechain=False)

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
