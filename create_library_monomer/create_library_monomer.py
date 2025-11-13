import glob
import os
import datetime
import shutil
import subprocess
import re
from topology.readmol.readpdbformat import ReadPdbFormat
from create_library_monomer_lib.parse_arguments import parse_arguments, print_header, command_info
from create_library_monomer_lib.extract_residues_from_itp import extract_residues_from_itp
from create_library_monomer_lib.get_all_ff_for_residues import get_all_ff_for_residues
from create_library_monomer_lib.extract_gromacs_parameters import extract_gromacs_parameters
from replicate_polymer_lib.logger import init_logger
from replicate_polymer_lib.check_connect_pdb import check_conect_pdb, check_resname_pdb, check_and_remove_ter_labels
from replicate_polymer_lib.remove_hydrogens import remove_hydrogens
from replicate_polymer_lib.typing_molecule import typing_molecule, write_gromacs_itps_xml, \
    write_gromacs_itps_itp, check_pdb_itp_consistency
from replicate_polymer_lib.insert_nonopenmm_ff_terms import insert_nonopenmm_ff_terms
from create_library_monomer_lib.gromacs_template import gmx_minim_template, gmx_eq1_nvt_template, \
                                                        gmx_eq2_npt_template, gmx_prod_npt_template


# =============================================================================
def function_pdb_type_ff_withxml(ipdb, opts, logger=None):

    # Read the original pdb and check PDB Connect section in the PDB file
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\tCheck resname and connect in the input file. ({})... ({})".format(ipdb, now)
    m += "\n\t\t{}".format(len(m)*"*")
    print(m) if logger is None else logger.info(m)
    filenamepdb = check_resname_pdb(ipdb)
    filenamepdb, bond_list_more100k = check_conect_pdb(filenamepdb, logger=logger)
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\tEnd Check resname and connect in the input file. ({})... ({})".format(ipdb, now)
    m += "\n\t\t{}".format(len(m)*"*")
    print(m) if logger is None else logger.info(m)

    # Remove hydrogen atoms
    if opts.hydrogens:
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\n\t\tRemoving hydrogens in the pdb file of the seed molecules.({})... ({})".format(filenamepdb, now)
        print(m) if logger is None else logger.info(m)
        filenamepdb_new, filenamegro_new = remove_hydrogens(filenamepdb, bond_list_more100k)
        filenamepdb = filenamepdb_new

    filenamepdb = check_and_remove_ter_labels(filenamepdb, logger=logger)

    # Build a ReadPdbFormat object for the molecular topology
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m1 = "\n\t\tTyping molecule ({})\n".format(now)
    m = "\t\t" + len(m1) * "*" + "\n"
    m += "\t\tLoading seed structure from {} ({})".format(filenamepdb, now)
    print(m1+m) if logger is None else logger.info(m1+m)
    untyped_mol = ReadPdbFormat(filenamepdb, logger=logger)

    # Apply force field to the seed molecules
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\tApplying force field to the seed molecules {} ({})".format(filenamepdb, now)
    print(m) if logger is None else logger.info(m)
    basepath = os.path.splitext(os.path.split(filenamepdb)[-1])[0]
    non_openmm_potentials_terms, lj14scale, coulomb14scale, ff = \
        typing_molecule(basepath+".top", untyped_mol, opts.forcefield, logger=logger,
                        type_kwargs={'verbose': opts.verbose,
                                     'assert_dihedral_params': True})

    # Insert forcefield terms which cannot be captured with openmm (i.e. Toxwaerd)
    insert_nonopenmm_ff_terms(filenamepdb, non_openmm_potentials_terms)

    isitp = True
    if isitp:
        base, _ = os.path.splitext(filenamepdb)
        base = os.path.split(base)[-1]
        write_gromacs_itps_xml(base+".top", ff, logger)

    # Generate itp and ff files for each residue
    itp_list = sorted(glob.glob(os.path.join("./forcefield", "*.itp")))
    defaults_dict, bondtype_dict, angletype_dict, dihedraltype_dict, dihedralfunc_dict, impropertype_avg_dict, \
        improper_idx_list, atomtypes = extract_gromacs_parameters(itp_list, [], [], wk=opts.pattern, log=logger)

    for itpfile in itp_list:
        residues_dict, itp_graph = extract_residues_from_itp(itpfile, wk=opts.pattern, debug=opts.debug)
        if itp_graph is not None:
            get_all_ff_for_residues(itp_graph, atomtypes,
                                    bondtype_dict, angletype_dict,
                                    dihedraltype_dict, dihedralfunc_dict, impropertype_avg_dict, improper_idx_list,
                                    defaults_dict, lj14scale,
                                    wk=opts.pattern,
                                    debug=opts.debug, log=logger)


# =============================================================================
def function_pdb_type_ff_withitp(ipdb, opts, list_ff_files, lj14scale, logger=None):

    # Read the original pdb and check PDB Connect section in the PDB file
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\tCheck resname and connect in the input file. ({})... ({})".format(ipdb, now)
    m += "\n\t\t{}".format(len(m)*"*")
    print(m) if logger is None else logger.info(m)
    filenamepdb = check_resname_pdb(ipdb)
    filenamepdb, bond_list_more100k = check_conect_pdb(filenamepdb, logger=logger)
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\tEnd Check resname and connect in the input file. ({})... ({})".format(ipdb, now)
    m += "\n\t\t{}".format(len(m)*"*")
    print(m) if logger is None else logger.info(m)

    filenamepdb = check_and_remove_ter_labels(filenamepdb, logger=logger)

    # Build a ReadPdbFormat object for the molecular topology
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m1 = "\n\t\tTyping molecule ({})\n".format(now)
    m = "\t\t" + len(m1) * "*" + "\n"
    m += "\t\tLoading seed structure from {} ({})".format(filenamepdb, now)
    print(m1+m) if logger is None else logger.info(m1+m)
    untyped_mol = ReadPdbFormat(filenamepdb, logger=logger)

    top_indices = [index for index, file in enumerate(list_ff_files) if 'top' in file]
    if len(top_indices) > 1:
        m = "\t\t ERROR: More than one top files in the list of forcefield files.\n"
        for idx in top_indices:
            m += "\t\t        {}\n".format(list_ff_files[idx])
        print(m) if logger is None else logger.info(m)
        exit()
    elif len(top_indices) == 0:
        # There is not top files. The first itp containing [ atoms ] and [ bonds ]
        labels = ['atoms', 'bonds']
        patterns = [re.compile(r'\[\s*' + label + r'\s*\]') for label in labels]
        for itp in list_ff_files:
            with open(itp, 'r') as fitp:
                lines = fitp.read()
                if all(pattern.search(lines) for pattern in patterns):
                    top_ffname = itp
                    break
    else:
        top_ffname = list_ff_files[top_indices[0]]

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m1 = "\n\t\tCheck for consistency between {} and {} ({})\n". \
        format(filenamepdb, top_ffname, now)
    m = "\t\t" + len(m1) * "*" + ""
    print(m1 + m) if logger is None else logger.info(m1 + m)
    # Check that the graph from pdb and itp are the same.
    # The order of the atoms and connectivity must be the same
    check, all_ffnames, itp_ffname = check_pdb_itp_consistency(untyped_mol, filenamepdb, top_ffname, logger)
    if not check:
        m = "\n\t\t ERROR. The PDB file is not consistent with the TOP file.\n"
        m += "\t\t       TOP file: {}\n".format(top_ffname)
        m += "\t\t       PDB file: {}".format(filenamepdb)
        print(m) if logger is None else logger.error(m)
        exit()
    if opts.noavgcharges:
        avg_charge = False
    else:
        avg_charge = True
    if opts.zerochargesres:
        zero_charge_res = True
    else:
        zero_charge_res = False

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m1 = "\n\t\tPreparing force field files from itp and pdb files ({})\n".format(now)
    m = "\t\t" + len(m1) * "*" + ""
    print(m1 + m) if logger is None else logger.info(m1 + m)
    is_write_pairs = write_gromacs_itps_itp(all_ffnames, itp_ffname, residuefile=opts.residues,
                                            chargefile=opts.charges,
                                            is_avg_charge=avg_charge,
                                            is_zero_charge_res=zero_charge_res, logger=logger)

    # Generate itp and ff files for each residue
    itp_list = sorted(glob.glob(os.path.join("./forcefield", "*.itp")))
    defaults_dict, bondtype_dict, angletype_dict, dihedraltype_dict, impropertype_dict, \
        improper_idx_list, atomtypes, dihedralfunc_dict = \
        extract_gromacs_parameters(itp_list, [], [], wk=opts.pattern, log=logger)

    for itpfile in itp_list:
        residues_dict, itp_graph = extract_residues_from_itp(itpfile, wk=opts.pattern, debug=opts.debug)
        if itp_graph is not None:
            get_all_ff_for_residues(itp_graph, atomtypes,
                                    bondtype_dict, angletype_dict,
                                    dihedraltype_dict, dihedralfunc_dict, impropertype_dict, improper_idx_list,
                                    defaults_dict, lj14scale, is_write_pairs,
                                    wk=opts.pattern,
                                    debug=opts.debug, log=logger)


# =============================================================================
def run_polyply_gen_itp(opts, logger):

    # Run the polyply program ==========================================
    sequences = opts.sequences
    lib_path = os.path.abspath(opts.pattern)
    namesystems = list()
    if opts.namesystems is None:
        for ipdb in opts.pdbfile:
            basename = os.path.splitext(os.path.split(ipdb)[-1])[0]
            namesystems.append(basename)
    else:
        namesystems = opts.namesystems

    for idx, ipdb in enumerate(opts.pdbfile):
        bash_command = "polyply gen_itp -lib {0:s} -seq  {1:s} -o {2:s}.itp -name {2:s}".\
            format(lib_path, sequences[idx], namesystems[idx])
        process = subprocess.Popen(bash_command.split(),
                                   stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\n\t\tOutput from polyply gen_itp. {} ({})".format(ipdb, now)
        m += "\n\t\t{}\n".format(len(m) * "*")
        m += output.decode()
        m += error.decode()
        print(m) if logger is None else logger.info(m)

    return namesystems


# =============================================================================
def run_polyply_gen_coords(namesystems, opts, logger):

    # Write general topology file
    # Write gromacs top for polyply
    filenametop_full = opts.pattern + ".top"
    with open(filenametop_full, 'w') as ftop:
        line = '#include "forcefield/forcefield.itp"\n'
        for ifile in namesystems:
            line += '#include "{}.itp"\n'.format(ifile)
        line += "\n"
        line += '[ system ]\n'
        line += '; name\n'
        line += 'Polymer\n\n'
        line += '[ molecules ]\n'
        line += ';name number\n'
        for idx, ifile in enumerate(namesystems):
            line += '{} {}\n'.format(ifile, opts.nmols[idx])

        ftop.writelines(line)

    # # gen_coords -p polymer.top -dens 650 -v -mir 1000 -mf 5000 -nr 10
    density, mir, mf, nr = opts.gen_coords_opts[0:]
    bash_command = "polyply gen_coords -p {0:s} -dens {1:f} -v -mir {2:d} -mf {3:d} -nr {4:d}". \
        format(filenametop_full, float(density), int(mir), int(mf), int(nr))
    process = subprocess.Popen(bash_command.split(),
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\tOutput from polyply gen_coord. {} ({})".format(filenametop_full, now)
    m += "\n\t\t{}\n".format(len(m) * "*")
    m += output.decode()
    m += error.decode()
    print(m) if logger is None else logger.info(m)


# =============================================================================
def run_gmx_minimization(opts, logger=None):

    # Create directory DIR_TARGET --> 01-MD_MIN
    # Check if the directory exists before attempting to delete it ======================
    dir_target = "./01-MD_MIN"
    if os.path.exists(dir_target):
        shutil.rmtree(dir_target)
    os.mkdir(dir_target)
    os.chdir(dir_target)

    # Generate minim.mdp
    gmx_minim_template(opts.gmxpath, opts.pattern + ".top", nsteps_min=50000)

    # RUN gromacs grompp
    grompp_exe = "{} grompp".format(opts.gmxpath)
    bash_command = "{} -f minim.mdp -c ../coords.gro -p ../{}". \
        format(grompp_exe, opts.pattern+".top")
    process = subprocess.Popen(bash_command.split(),
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\tOutput from grompp Gromacs.({})".format(now)
    m += "\n\t\t{}\n".format(len(m) * "*")
    m += output.decode()
    m += error.decode()
    print(m) if logger is None else logger.info(m)

    # RUN gromacs mdrun
    mdrun_exe = "{} mdrun".format(opts.gmxpath)
    bash_command = "{}". \
        format(mdrun_exe)
    process = subprocess.Popen(bash_command.split(),
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\tOutput from mdrun Gromacs.({})".format(now)
    m += "\n\t\t{}\n".format(len(m) * "*")
    m += output.decode()
    m += error.decode()
    print(m) if logger is None else logger.info(m)

    # RUN gromacs trjconv
    trjconv_exe = "echo 0|{} trjconv".format(opts.gmxpath)
    bash_command = "{} -f traj.trr -o traj_unwrap.trr -pbc whole". \
        format(trjconv_exe)
    process = subprocess.Popen(bash_command, shell=True,
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\tOutput from trjconv Gromacs.({})".format(now)
    m += "\n\t\t{}\n".format(len(m) * "*")
    m += output.decode()
    m += error.decode()
    print(m) if logger is None else logger.info(m)

    os.chdir("../")


# =============================================================================
def main_app(version):

    # Init logger =======================================================================
    logger = init_logger("Output", fileoutput="Info.log", append=False, inscreen=True)
    print_header(version, logger)

    # Parse arguments in the command line ===============================================
    opts = parse_arguments()

    # Info ==============================================================================
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    starttime = datetime.datetime.now()
    m = "\t\tStart Job at {} ============".format(now)
    print(m) if logger is None else logger.info(m)

    # Print command info ================================================================
    command_info(opts, logger=logger)

    # Check files in the force field directory ==========================================
    if os.path.isfile(opts.forcefield):
        ff_extension = os.path.splitext(os.path.split(opts.forcefield)[-1])[-1]
        list_files = [opts.forcefield]
    else:
        ff_dir = opts.forcefield
        list_files_top = glob.glob(os.path.join(ff_dir, '*.top'))
        list_files_itp = glob.glob(os.path.join(ff_dir, '*.itp'))
        list_files = list_files_top + list_files_itp
        if len(list_files_top) != 0:
            ff_extension = ".top"
        else:
            ff_extension = ".itp"

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\tInitial FF Files. ({})\n".format(now)
    for item in list_files:
        m += "\t\t\t{}\n".format(item)
    print(m) if logger is None else logger.info(m)

    # Delete previous generated files ===================================================
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    if opts.nodeleteff:
        m = "\t\tFiles with extension ff for PolyPly have been deleted. ({})".format(now)
        for ifile in glob.glob("*.ff"):
            os.remove(ifile)
    else:
        m = "\t\tFiles with extension ff for PolyPly have not been deleted. ({})".format(now)
    print(m) if logger is None else logger.info(m)

    # Create directory DIR_TARGET
    # Check if the directory exists before attempting to delete it ======================
    dir_target = opts.pattern
    if opts.pattern is not None and os.path.exists(dir_target):
        shutil.rmtree(dir_target)
    os.mkdir(dir_target)
    # Loop over the pdbs
    for ipdb in opts.pdbfile:
        if ff_extension == ".xml":
            function_pdb_type_ff_withxml(ipdb, opts, logger=logger)
        elif ff_extension == ".itp" or ff_extension == ".top":
            lj14scale = 0.5
            function_pdb_type_ff_withitp(ipdb, opts, list_files, lj14scale, logger=logger)

    # Run polyply gen_itp ===============================================================
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    if opts.sequences is not None:
        m = "\n\t\tStart polyply to generate ITP files. ({})\n".format(now)
        print(m) if logger is None else logger.info(m)
        namesystems = run_polyply_gen_itp(opts, logger)
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\n\t\tEnd polyply to generate ITP files. ({})\n".format(now)
        print(m) if logger is None else logger.info(m)
    else:
        namesystems = []

    # Run polyply gen_coords ============================================================
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    if opts.sequences is not None and opts.nmols is not None:
        m = "\n\t\tStart polyply to generate Coordinate files. ({})\n".format(now)
        print(m) if logger is None else logger.info(m)
        run_polyply_gen_coords(namesystems, opts, logger)
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\n\t\tEnd polyply to generate Coordinate files. ({})\n".format(now)
        print(m) if logger is None else logger.info(m)

    # Run gromacs minimization ===========================================================
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    if opts.gmxpath is not None:
        if os.path.isfile(opts.gmxpath):
            m = "\n\t\tStart gromacs minimization. ({})\n".format(now)
            print(m) if logger is None else logger.info(m)
            run_gmx_minimization(opts, logger=logger)
            now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
            m = "\n\t\tEnd gromacs minimization. ({})\n".format(now)
            print(m) if logger is None else logger.info(m)
            # Genarate templates for equilibration
            m = "\n\t\tGenerate mdp and run templates for GROMACS. ({})\n".format(now)
            print(m) if logger is None else logger.info(m)
            gmx_eq1_nvt_template(opts.gmxpath, opts.pattern + ".top")
            gmx_eq2_npt_template(opts.gmxpath, opts.pattern + ".top")
            gmx_prod_npt_template(opts.gmxpath, opts.pattern + ".top")
            m = "\n\t\tEnd Generate mdp and run templates for GROMACS. ({})\n".format(now)
            print(m) if logger is None else logger.info(m)
        else:
            m = "\n\t\t Path to GROMACS does not exist.\n"
            m += "\n\t\t {}\n".format(opts.gmxpath)
            print(m) if logger is None else logger.warn(m)

    # Move the force field files
    dir_target = "./00-FORCE_FIELD_POLYPLY"
    if os.path.exists(dir_target):
        shutil.rmtree(dir_target)
    os.mkdir(dir_target)
    try:
        shutil.move("./"+opts.pattern, dir_target)
        shutil.move("forcefield", dir_target)
        shutil.move(opts.pattern+".top", dir_target)
        for iname in namesystems:
            shutil.move(iname + ".itp", dir_target)
    except FileNotFoundError:
        pass

    # End of Job ========================================================================
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
