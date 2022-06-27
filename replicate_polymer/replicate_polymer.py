import os
import sys
import datetime
from  topology.readmol.readpdbformat import ReadPdbFormat
from replicate_polymer_lib.parse_arguments import parse_arguments, print_header, command_info
from replicate_polymer_lib.logger import init_logger
from replicate_polymer_lib.check_connect_pdb import check_conect_pdb, check_resname_pdb
from replicate_polymer_lib.remove_hydrogens import remove_hydrogens
from replicate_polymer_lib.typing_molecule import typing_molecule
from replicate_polymer_lib.replicate_pdb import replicate_pdb
from replicate_polymer_lib.prepare_lammps import prepare_lammps, clean_lammps
from replicate_polymer_lib.exclusions import exclusion_gromacs, insert_exclusions
from replicate_polymer_lib.insert_impropers_top_gromacs import insert_impropers
from replicate_polymer_lib.insert_nonopenmm_ff_terms import insert_nonopenmm_ff_terms

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
    print(m) if logger is None else logger.info(m)
    filenamepdb = check_resname_pdb(opts.pdbfile)
    filenamepdb = check_conect_pdb(filenamepdb)

    # Remove hydrogen atoms
    if opts.hydrogens:
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\n\t\tRemoving hydrogens in the pdb file of the seed molecles.({})... ({})".format(filenamepdb, now)
        print(m) if logger is None else logger.info(m)
        filenamepdb_new, filenamegro_new = remove_hydrogens(filenamepdb)
        filenamepdb = filenamepdb_new

    # Build a ReadPdbFormat object for the molecular topology
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m1 = "\t\tReplicating and typing molecule ({})\n".format(now)
    m = "\t\t" + len(m1) * "*" + "\n"
    m += "\t\tLoading seed structure from {} ({})".format(filenamepdb, now)
    print(m1+m) if logger is None else logger.info(m1+m)
    untyped_mol = ReadPdbFormat(filenamepdb)

    # Apply force field to the seed molecules
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\tApplying force field to the seed molecules {} ({})".format(filenamepdb, now)
    print(m) if logger is None else logger.info(m)
    basepath = os.path.splitext(os.path.split(filenamepdb)[-1])[0]
    non_openmm_potentials_terms = typing_molecule(basepath+".top", untyped_mol, opts.xmlfile, logger=logger,
                                  type_kwargs={'verbose': opts.verbose, 'assert_dihedral_params': True})

    # Replicate the system (image)
    if opts.images is not None:
        m = "\t\tReplicating molecules {} ({})".format(filenamepdb, now)
        print(m) if logger is None else logger.info(m)
        _, _, residues_map_atom = replicate_pdb(filenamepdb, opts.images, boxlength=opts.boxlength,
                                                            boxangles=opts.boxangle, logger_log=logger)

    # Insert improper angles in the top file
    insert_impropers(opts, filenamepdb, logger=logger)

    # Insert exclusions in the top file
    insert_exclusions(opts, filenamepdb, residues_map_atom, logger=logger)

    # Insert forcefield terms which cannot be captured with openmm (i.e. Toxwaerd)
    insert_nonopenmm_ff_terms(filenamepdb, non_openmm_potentials_terms)

    # If engine is Lammps also lammps files are generated using (intermol, https://github.com/shirtsgroup/InterMol)
    if opts.engine == "lammps":
        base = os.path.splitext(os.path.split(filenamepdb)[-1])[0]
        basefile_gro = basepath+'_replicate.gro'
        basefile_top = basepath+'_replicate.top'
        dirlocal = "./"
        ext = ".gro"
        filenamegro_new = os.path.join(dirlocal, base + "_replicate" + ext)
        input_filename, data_filename = prepare_lammps([basefile_gro, basefile_top], title=filenamegro_new, logger=logger)
        clean_lammps(input_filename, data_filename, basefile_gro, basefile_top, opts.xmlfile)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    endtime = datetime.datetime.now()
    m = "\t\tJob  Done at {} ============\n".format(now)
    print(m) if logger is None else logger.info(m)
    m = "\t\tTotal time: {0:.2f} seconds".format((endtime-starttime).total_seconds())
    print(m) if logger is None else logger.info(m)



# =============================================================================
if __name__ == "__main__":

    import version
    main_app(version.__version__)
