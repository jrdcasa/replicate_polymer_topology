__version__ = "3.0"

from typing import Optional, List, TextIO


def main(argv: Optional[List[str]] = None, stream: Optional[TextIO] = None) -> int:

    from create_library_monomer.create_library_monomer import main_app

    #with contextlib.ExitStack() as ctx:
    return main_app(__version__)
