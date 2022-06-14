import os
import datetime
from replicate_polymer_lib.setup_chiral_impropers import setup_chiral_impropers


def insert_impropers(opts, fnamepdb, logger=None):

    base, _ = os.path.splitext(fnamepdb)
    base = os.path.split(base)[-1]
    dirlocal = "./"

    # Find improper angles
    if opts.xmlfile.find("trappe") != -1:
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\tInserting impropers in the topology file ({})\n".format(now)
        if opts.impropers.title() == "Guess":
            m += "\t\t\tGuessing impropers ..."
            nimpropers, imp_lines = setup_chiral_impropers(fnamepdb, improper_filename=None)

        else:
            nimpropers, imp_lines = setup_chiral_impropers(fnamepdb, improper_filename=opts.impropers)
            m += "\t\t\tInserting impropers from {} ...".format(opts.impropers)
        print(m) if logger is None else logger.info(m)
    else:
        imp_lines = ""
        nimpropers = 0

    if nimpropers != 0:
        ext = ".top"
        fname = os.path.join(dirlocal, base + "_replicate" + ext)
        with open(fname, "r") as f:
            contents = f.readlines()

        idx_to_insert = contents.index("[ dihedrals ]\n")

        contents.insert(idx_to_insert + 1, imp_lines)

        with open(fname, "w") as f:
            contents = "".join(contents)
            f.write(contents)
