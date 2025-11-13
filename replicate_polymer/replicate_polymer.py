import os
import sys
import datetime
from topology.readmol.readpdbformat import ReadPdbFormat
from replicate_polymer_lib.parse_arguments import parse_arguments, print_header, command_info
from replicate_polymer_lib.logger import init_logger
from replicate_polymer_lib.check_connect_pdb import check_conect_pdb, check_resname_pdb, check_and_remove_ter_labels
from replicate_polymer_lib.remove_hydrogens import remove_hydrogens
from replicate_polymer_lib.typing_molecule import typing_molecule
from replicate_polymer_lib.replicate_pdb import replicate_pdb
from replicate_polymer_lib.prepare_lammps import clean_lammps
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
    m += "\n\t\t{}".format(len(m)*"*")
    print(m) if logger is None else logger.info(m)
    filenamepdb = check_resname_pdb(opts.pdbfile)
    filenamepdb, bond_list_more100K = check_conect_pdb(filenamepdb, logger=logger)
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\tEnd Check resname and connect in the input file. ({})... ({})".format(opts.pdbfile, now)
    m += "\n\t\t{}".format(len(m)*"*")
    print(m) if logger is None else logger.info(m)

    # Remove hydrogen atoms
    if opts.hydrogens:
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\n\t\tRemoving hydrogens in the pdb file of the seed molecules.({})... ({})".format(filenamepdb, now)
        print(m) if logger is None else logger.info(m)
        filenamepdb_new, filenamegro_new = remove_hydrogens(filenamepdb, bond_list_more100K)
        filenamepdb = filenamepdb_new

    filenamepdb = check_and_remove_ter_labels(filenamepdb, logger=logger)

    # Build a ReadPdbFormat object for the molecular topology
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m1 = "\n\t\tReplicating and typing molecule ({})\n".format(now)
    m = "\t\t" + len(m1) * "*" + "\n"
    m += "\t\tLoading seed structure from {} ({})".format(filenamepdb, now)
    print(m1+m) if logger is None else logger.info(m1+m)
    untyped_mol = ReadPdbFormat(filenamepdb, logger=logger)

    # Apply force field to the seed molecules
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\tApplying force field to the seed molecules {} ({})".format(filenamepdb, now)
    print(m) if logger is None else logger.info(m)
    basepath = os.path.splitext(os.path.split(filenamepdb)[-1])[0]
    non_openmm_potentials_terms, lj14scale, coulomb14scale, _ = typing_molecule(basepath+".top", untyped_mol,
                                  opts.xmlfile, logger=logger,
                                  type_kwargs={'verbose': opts.verbose, 'assert_dihedral_params': True})

    # Replicate the system (image)
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    if opts.images is not None:
        m = "\t\tReplicating molecules {} ({})".format(filenamepdb, now)
        print(m) if logger is None else logger.info(m)
        _, _, residues_map_atom, mols_replicate = replicate_pdb(filenamepdb, opts.images, boxlength=opts.boxlength,
                                                                boxangles=opts.boxangle, logger_log=logger)

    # Insert forcefield terms which cannot be captured with openmm (i.e. Toxwaerd)
    base, _ = os.path.splitext(filenamepdb)
    base = os.path.split(base)[-1]
    insert_nonopenmm_ff_terms(base+"_replicate.pdb", non_openmm_potentials_terms)

    # Insert improper angles in the top file
    if opts.impropers.upper() != 'NO':
        insert_impropers(opts, base+"_replicate.pdb", logger=logger)
    elif opts.hydrogens and opts.impropers.upper():
        m = "\t\tWARNING: Using a united-atom model without guess or define improper angles. " \
            "\n\t\tIf this is correct ignore this warning."
        print(m) if logger is None else logger.info(m)


    # Insert exclusions in the top file
    insert_exclusions(opts, filenamepdb+"_replicate", residues_map_atom, logger=logger)

    # If engine is Lammps also lammps files are generated
    if opts.engine == "lammps" or opts.engine == "lammps_mc" :
        # Apply force field to the seed molecules
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\tCreating LAMMPS input files ({})".format(now)
        print(m) if logger is None else logger.info(m)
        base = os.path.splitext(os.path.split(filenamepdb)[-1])[0]
        basefile_gro = basepath+'_replicate.gro'
        basefile_top = basepath+'_replicate.top'
        dirlocal = "./"
        ext = ".gro"
        filenamegro_new = os.path.join(dirlocal, base + "_replicate" + ext)
        if opts.engine == "lammps_mc":
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
