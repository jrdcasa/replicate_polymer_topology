__version__ = "3.0"

from typing import Optional, List, TextIO


def main(argv: Optional[List[str]] = None, stream: Optional[TextIO] = None) -> int:

    from remove_hydrogens.remove_hydrogens import main_app

    #with contextlib.ExitStack() as ctx:
    return main_app(__version__)