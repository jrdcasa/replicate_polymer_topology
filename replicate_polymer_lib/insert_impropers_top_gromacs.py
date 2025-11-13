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
            m += "\t\t\tGuessing impropers ... (it is better use a improper list)"
            nimpropers_list, imp_lines_list = setup_chiral_impropers(fnamepdb, improper_filename=None)
        else:
            nimpropers_list, imp_lines_list = setup_chiral_impropers(fnamepdb, improper_filename=opts.impropers)
            m += "\t\t\tInserting impropers from {} ...".format(opts.impropers)
        print(m) if logger is None else logger.info(m)
    elif opts.impropers is not None:
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\tInserting impropers in the topology file ({})\n".format(now)
        nimpropers_list, imp_lines_list = setup_chiral_impropers(fnamepdb, improper_filename=opts.impropers,
                                                                 logging=logger)
        m += "\t\t\tInserting impropers from {} ...".format(opts.impropers)
        print(m) if logger is None else logger.info(m)
    else:
        imp_lines = [""]
        nimpropers_list = [0]

    # File with improper data
    ext = ".top"
    fname = os.path.join(dirlocal, base + ext)
    with open(fname, "r") as f:
        contents = f.readlines()

        # Molecule Type and atoms_labels
        idx_moleculetype_labels = []
        idx_atoms_labels = []
        idx_insert_return_atend = []
        idx_dihedral_insert = []
        for idx in range(len(contents)):
            if contents[idx] == "[ moleculetype ]\n":
                idx_moleculetype_labels.append(idx)
            if contents[idx] == "[ atoms ]\n":
                idx_atoms_labels.append(idx)
        idx_moleculetype_labels.append(len(contents))

        for idx in range(len(idx_moleculetype_labels)-1):
            istart = idx_moleculetype_labels[idx]
            iend = idx_moleculetype_labels[idx + 1]
            # Find the position to insert the new dihedrals
            pos = -1
            # Insert after [ dihedrals ] labels, if exists
            j = istart
            for iline in contents[istart:iend]:
                if iline == "[ dihedrals ]\n":
                    pos = j
                    idx_insert_return_atend.append(False)
                    break
                else:
                    j += 1
            # Insert after [ angles ] labels.
            if pos == -1:
                j = istart
                for iline in contents[istart:iend]:
                    if iline == "[ angles ]\n":
                        pos = j
                        for jline in contents[pos:iend]:
                            if len(jline) < 2:
                                pos = j + 1
                                contents.insert(pos, "[ dihedrals ]\n")
                                idx_insert_return_atend.append(True)
                                pos = j + 2
                                break
                            else:
                                j += 1
                    else:
                        j += 1

            idx_dihedral_insert.append(pos+1)

        for i in range(len(idx_dihedral_insert)):
            if nimpropers_list[i] != 0:
                idx_to_insert = idx_dihedral_insert[i]
                contents.insert(idx_to_insert+1, imp_lines_list[i])

        with open(fname, "w") as f:
            contents = "".join(contents)
            f.write(contents)
