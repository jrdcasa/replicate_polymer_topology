__version__ = "3.0"

import contextlib
from typing import Optional, List, TextIO


def main(argv: Optional[List[str]] = None, stream: Optional[TextIO] = None) -> int:

    from create_lammps_from_gromacs.create_lammps_from_gromacs import main_app

    #with contextlib.ExitStack() as ctx:
    return main_app(__version__)
