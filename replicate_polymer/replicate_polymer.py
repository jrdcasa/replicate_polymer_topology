import os
import sys
import datetime
from  topology.readmol.readpdbformat import ReadPdbFormat
from replicate_polymer_lib.parse_arguments import parse_arguments, print_header
from replicate_polymer_lib.logger import init_logger
from replicate_polymer_lib.check_connect_pdb import check_conect_pdb
from replicate_polymer_lib.remove_hydrogens import remove_hydrogens
from replicate_polymer_lib.setup_chiral_impropers import setup_chiral_impropers
from replicate_polymer_lib.typing_molecule import typing_molecule
from replicate_polymer_lib.replicate_pdb import replicate_pdb
from replicate_polymer_lib.prepare_lammps import prepare_lammps
from replicate_polymer_lib.exclusions import exclusion_gromacs

# =============================================================================
def main_app(version):

    # Parse arguments in the command line
    opts = parse_arguments()

    # Init logger
    logger = init_logger("Output", fileoutput="Info.log", append=False, inscreen=True)
    print_header(version, logger)
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

    m = "\t\tStart Job at {} ============".format(now)
    print(m) if logger is None else logger.info(m)

    m1 = ""
    for item in sys.argv[1:]:
        m1 += " {}".format(item)
    m = "\n\t\tCommand line: \n"
    m += "\t\t\tpython {}".format(os.path.split(sys.argv[0])[1])
    m += m1+"\n"
    m += "\t\t\t         or\n"
    m += "\t\t\treplicate_polymer".format(os.path.split(sys.argv[0])[1])
    m += m1+"\n"
    print(m) if logger is None else logger.info(m)

    # Read the original pdb and check PDB Connect section in the PDB file
    filenamepdb = check_conect_pdb(opts.pdbfile)
    if filenamepdb is None:
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\n\t\tRemoving hydrogens in the pdb file of the seed molecles.({})... ({})".format(filenamepdb, now)
        print(m) if logger is None else logger.info(m)
        exit()

    # Remove hydrogen atoms
    if opts.hydrogens:
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\n\t\tRemoving hydrogens in the pdb file of the seed molecles.({})... ({})".format(filenamepdb, now)
        print(m) if logger is None else logger.info(m)
        filenamepdb_new, filenamegro_new = remove_hydrogens(filenamepdb)
        filenamepdb = filenamepdb_new

    # IMPROPERS for Trappe-UA force field and chiral atoms
    if opts.xmlfile.find("trappe") != -1:
        if opts.impropers.title() == "Guess":
            nimpropers, imp_lines = setup_chiral_impropers(filenamepdb, improper_filename=None)
        else:
            nimpropers, imp_lines = setup_chiral_impropers(filenamepdb, improper_filename=opts.impropers)
    else:
        imp_lines = ""
        nimpropers = 0

    # Build the Compound object of the seed molecules from mbuild
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m1 = "\t\tReplicating and typing molecule ({})\n".format(now)
    m = "\t\t" + len(m1) * "*" + "\n"
    m += "\t\tLoading seed structure from {} ({})".format(filenamepdb, now)
    print(m1+m) if logger is None else logger.info(m1+m)
    untyped_mol = ReadPdbFormat(filenamepdb)

    # Apply force field to the seed molecule
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\tApplying force field to the seed molecules {} ({})".format(filenamepdb, now)
    print(m) if logger is None else logger.info(m)
    basepath = os.path.splitext(os.path.split(filenamepdb)[-1])[0]
    # Each chain (or molecule) is a residue
    typing_molecule(basepath+".top", untyped_mol, opts.xmlfile, logger=logger,
                    type_kwargs={'verbose': True, 'assert_dihedral_params': True})

    # Replicate the pdb chain
    if opts.images is not None:
        _, trj_replicate, residues_map_atom = replicate_pdb(filenamepdb, opts.images, boxlength=opts.boxlength,
                                                            boxangles=opts.boxangle, logger_log=logger)

    base, _ = os.path.splitext(filenamepdb)
    base = os.path.split(base)[-1]
    dirlocal = "./"
    ext = ".gro"
    filenamegro_new = os.path.join(dirlocal, base + "_replicate" + ext)

    # Insert improper angles in the top file
    ext = ".top"
    fname = os.path.join(dirlocal, base + "_replicate" + ext)
    with open(fname, "r") as f:
         contents = f.readlines()

    idx_to_insert = contents.index("[ dihedrals ]\n")

    contents.insert(idx_to_insert+1, imp_lines)

    with open(fname, "w") as f:
        contents = "".join(contents)
        f.write(contents)

    # Insert exclusions in the top file
    if opts.npairs is not None:
        ext = ".top"
        fname = os.path.join(dirlocal, base + "_replicate" + ext)
        with open(fname, "r") as f:
            contents = f.readlines()

        idx_to_insert = contents.index("[ system ]\n")

        lines_exclusions = exclusion_gromacs(opts.npairs, residues_map_atom)

        contents.insert(idx_to_insert-1, lines_exclusions)

        with open(fname, "w") as f:
            contents = "".join(contents)
            f.write(contents)

    # If engine is Lammps also lammps files are generated using (intermol, https://github.com/shirtsgroup/InterMol)
    if opts.engine == "lammps":
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m1 = "\t\tGenerating lmp and data files for {} ({})\n".format("LAMMPS", now)
        m = "\t\t" + len(m1) * "*" + "\n"
        print(m1+m) if logger is None else logger.info(m1+m)
        basepath = os.path.splitext(os.path.split(filenamegro_new)[-1])[0]
        basefile_gro = basepath+'.gro'
        basefile_top = basepath+'.top'
        prepare_lammps([basefile_gro, basefile_top], title=filenamegro_new)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\tJob  Done at {} ============\n".format(now)
    print(m) if logger is None else logger.info(m)


# =============================================================================
if __name__ == "__main__":

    import version
    main_app(version.__version__)
