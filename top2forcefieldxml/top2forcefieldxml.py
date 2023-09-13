from replicate_polymer_lib.logger import init_logger
from top2forcefieldxml_lib.parse_arguments import print_header


# =============================================================================
def main_app(version):

    # Init logger
    logger = init_logger("Output", fileoutput="Info_top.log", append=False, inscreen=True)
    print_header(version, logger)


# =============================================================================
if __name__ == "__main__":

    __version__ = "1.1"
    main_app(__version__)
    