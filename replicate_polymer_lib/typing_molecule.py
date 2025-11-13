import numpy as np
import datetime
import re
import os
import shutil
import statistics
from replicate_polymer_lib.forcefield.Forcefield import Forcefield
from create_library_monomer_lib.extract_itp_parameters import atomic_mass, atomic_number
from collections import defaultdict


# =================================================================================================
def typing_molecule(filename, molecule, forcefield_xml_files, overwrite=True,
                    type_kwargs=None, logger=None, **kwargs):
    """

    Args:
        filename: Name of the file to write the topology using a Parmed instance
        molecule: A ReadPDBFormat object from topology library
        forcefield_xml_files: List of XML files containing force field parameters
        overwrite : bool, optional, default=False
                    Overwrite if the filename already exists
        type_kwargs:  dict, optional, default=None
                    Keyword arguments to provide to `Forcefield.apply`.
        logger: Logger object

    Returns:

    """

    # Parse force field xml file
    ff = Forcefield(forcefield_files=forcefield_xml_files, debug=False, logger=logger)

    # Extra arguments
    if not type_kwargs:
        type_kwargs = {}

    # Non-periodic systems --> Create a box
    is_all_zero = np.all((molecule._unitcell == 0))
    if is_all_zero:
        xmin = min(molecule._universe.coord.positions[:, 0])
        ymin = min(molecule._universe.coord.positions[:, 1])
        zmin = min(molecule._universe.coord.positions[:, 2])
        xmax = max(molecule._universe.coord.positions[:, 0])
        ymax = max(molecule._universe.coord.positions[:, 1])
        zmax = max(molecule._universe.coord.positions[:, 2])
        molecule._unitcell = [[xmax - xmin + 5., 0.0, 0.0],
                              [0.0, ymax - ymin + 5., 0.0],
                              [0.0, 0.0, zmax - zmin + 5.]]
        molecule._boxlength = [xmax - xmin + 5., ymax - ymin + 5., zmax - zmin + 5.]
        molecule._boxangle = [90.0, 90.0, 90.0]

    pmd_structure, non_openmm_potentials_terms = ff.apply(molecule, logger=logger, **type_kwargs)

    total_charge = sum([atom.charge for atom in pmd_structure])
    if round(total_charge, 4) != 0.0:
        m1 = ("\t\tWARN: System is not charge neutral. Total charge is {}.\n".format(total_charge))
        m2 = "\t\tCheck charges!!!!!\n"
        ll = int(len(m1) / 1)
        m = "\t\t" + ll * "-" + "\n"
        print(m + m1 + m2 + m) if logger is None else logger.info(m + m1 + m2 + m)

    import sys
    sys.setrecursionlimit(30000)
    # Parmed sets gen_pairs to 'no' so that the pair-specific L-J parameters are
    # printed to the topology file (rather than being auto-created)
    pmd_structure.save(filename, overwrite=overwrite, **kwargs)

    # Recover lj14scale and coulomb14scale from xml file
    with open(forcefield_xml_files, 'r') as fxml:
        for line in fxml:
            if line.count("lj14scale") != 0 and line.count("coulomb14scale") == 0:
                # Using regular expression to extract values between quotes
                matches = re.findall(r'lj14scale="(.*?)"', line)
                lj14scale = float(matches[0])
                coulomb14scale = 1.0
            elif line.count("lj14scale") == 0 and line.count("coulomb14scale") != 0:
                # Using regular expression to extract values between quotes
                matches = re.findall(r'coulomb14scale="(.*?)"', line)
                lj14scale = 1.0
                coulomb14scale = float(matches[0])
            elif line.count("lj14scale") != 0 and line.count("coulomb14scale") != 0:
                # Using regular expression to extract values between quotes
                matches = re.findall(r'coulomb14scale="(.*?)"', line)
                coulomb14scale = float(matches[0])
                matches = re.findall(r'lj14scale="(.*?)"', line)
                lj14scale = float(matches[0])

    # if isitp:
    #     write_gromacs_itps_itp(filename, ff, logger)

    return non_openmm_potentials_terms, lj14scale, coulomb14scale, ff


# =================================================================================================
def write_gromacs_itps_xml(filename_top, ff, logger=None):
    # Create directory forcefield
    dir_target = "forcefield"
    if os.path.exists(dir_target):
        shutil.rmtree(dir_target)
    os.mkdir(dir_target)

    with open(filename_top, 'r') as ftop:
        lines = ftop.readlines()
        # Get each molecule
        pos_moleculetype = [index for index, iline in enumerate(lines) if 'moleculetype' in iline.lower()]
        ntypemols = len(pos_moleculetype)
        pos_system = [index for index, iline in enumerate(lines) if re.match(r"^\[(\s*)system(\s*)]", iline)]
        pos_moleculetype.append(pos_system[0])
        # Get defaults label position
        pos_defaults = [index for index, iline in enumerate(lines) if re.match(r"^\[(\s*)defaults(\s*)]", iline)]
        # Get atomtypes
        pos_atomtypes = [index for index, iline in enumerate(lines) if re.match(r"^\[(\s*)atomtypes(\s*)]", iline)]
        # Get name and molecules
        # pos_molecules = [index for index, iline in enumerate(lines) if re.match(r"^\[(\s*)molecules(\s*)]", iline)]

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m1 = "\n\t\tGenerating itps for molecules in directory forcefield ({})\n".format(now)
        m = "\t\t" + len(m1) * "*"
        print(m1 + m) if logger is None else logger.info(m1 + m)

        # Generate the itp files
        itp_lines = ""
        for imol in range(ntypemols):
            idx_start = pos_moleculetype[imol]
            idx_end = pos_moleculetype[imol + 1]
            # Get name of the molecules
            idx = idx_start
            while True:
                iline = lines[idx]
                if len(iline.replace(" ", "")) < 2:
                    break
                if not re.match("^;", iline.strip()):
                    name_mol = iline.split()[0]
                idx += 1
            itp_lines += '#include "{}.itp"\n'.format(name_mol)
            m1 = "\t\t\tGenerating ./forcefield/{}.itp".format(name_mol)
            print(m1) if logger is None else logger.info(m1)
            with open(os.path.join(dir_target, name_mol + ".itp"), 'w') as fitp:
                line_header = "; Topology file write with create_library_monomer toolkit\n"
                line_header += ";                Javier Ramos (IEM-CSIC)\n"
                fitp.writelines(line_header)
                for idx in range(idx_start, idx_end):
                    fitp.writelines(lines[idx])

        m1 = "\t\t\tGenerating ./forcefield/forcefield.itp and ffnonbonded.itp"
        print(m1) if logger is None else logger.info(m1)

        # Generate forcefield.itp and ffnonbonded.itp files
        with open(os.path.join(dir_target, "forcefield.itp"), 'w') as fff:
            idx = pos_defaults[0]
            while True:
                iline = lines[idx]
                if len(iline.replace(" ", "")) < 2:
                    break
                fff.writelines(iline)
                idx += 1
            fff.writelines("\n")
            fff.writelines('#include "ffnonbonded.itp"')

        # Generate forcefield.itp and ffnonbonded.itp files
        with open(os.path.join(dir_target, "ffnonbonded.itp"), 'w') as fff:
            idx = pos_atomtypes[0]
            while True:
                iline = lines[idx]
                if len(iline.replace(" ", "")) < 2:
                    break
                try:
                    label = iline.split()[0]
                    iline = iline.split()[0] + \
                        " {0:s} ".format(ff.atomTypeClasses[label]) + \
                        '  '.join(map(str, iline.split()[1:])) + "\n"
                except KeyError:
                    pass
                fff.writelines(iline)
                idx += 1


# =================================================================================================
def write_gromacs_itps_itp(all_ffnames, itp_ff, chargefile=None, residuefile=None,
                           is_avg_charge=True, is_zero_charge_res=False, logger=None):

    # =========================================================================
    def aux_find_end(lines_aux, start):

        kdx_list = list()
        if len(start) != 0:
            for kdx in start:
                while True:
                    try:
                        iline = lines_aux[kdx+1]
                    except IndexError:
                        break
                    if len(iline.replace(" ", "")) < 2:
                        break
                    if re.match(r'^\[', iline.strip()):
                        break
                    kdx += 1
                kdx_list.append(kdx)
        else:
            kdx_list.append(-1)

        return kdx_list

    # ==========================================================================
    def aux_ntokens(label_aux, pos_dict_start):

        ntoken_list = list()
        for ikey, item in pos_dict_start.items():
            if len(item) != 0:
                ntoken_list.append(ikey)
        if len(ntoken_list) > 1:
            m = "\n\t\t ERROR. {} label appears more than one time.\n".format(label_aux)
            for i in ntoken_list:
                m += "\t\t       {} in {}\n".format(label_aux, i)
            print(m) if logger is None else logger.error(m)
            exit()

        return ntoken_list

    # ==========================================================================
    # Create directory forcefield
    dir_target = "forcefield"
    if os.path.exists(dir_target):
        shutil.rmtree(dir_target)
    os.mkdir(dir_target)

    pos_include = defaultdict(list)
    pos_system_start = defaultdict(list)
    pos_molecules_start = defaultdict(list)
    pos_moleculetype_start = defaultdict(list)
    pos_defaults_start = defaultdict(list)
    pos_atomtypes_start = defaultdict(list)
    pos_bondtypes_start = defaultdict(list)
    pos_pairtypes_start = defaultdict(list)
    pos_angletypes_start = defaultdict(list)
    pos_dihedraltypes_start = defaultdict(list)
    pos_atoms_start = defaultdict(list)
    pos_bonds_start = defaultdict(list)
    pos_pairs_start = defaultdict(list)
    pos_angles_start = defaultdict(list)
    pos_dihedrals_start = defaultdict(list)

    pos_include_end = defaultdict(list)
    pos_system_end = defaultdict(list)
    pos_molecules_end = defaultdict(list)
    pos_moleculetype_end = defaultdict(list)
    pos_defaults_end = defaultdict(list)
    pos_atomtypes_end = defaultdict(list)
    pos_bondtypes_end = defaultdict(list)
    pos_pairtypes_end = defaultdict(list)
    pos_angletypes_end = defaultdict(list)
    pos_dihedraltypes_end = defaultdict(list)
    pos_atoms_end = defaultdict(list)
    pos_bonds_end = defaultdict(list)
    pos_pairs_end = defaultdict(list)
    pos_angles_end = defaultdict(list)
    pos_dihedrals_end = defaultdict(list)

    # Read the start and end positions of each file containing force field information
    for ifile in all_ffnames:
        # Read the itp and write it
        with open(ifile, 'r') as fitp_ff:

            # Read the itp and write it
            lines = fitp_ff.readlines()
            ifile_name = ifile
            # Include
            pos_include[ifile_name] = [index for index, iline in enumerate(lines)
                                       if '#include' in iline.lower()]
            # [ system ]
            pos_system_start[ifile_name] = [index for index, iline in enumerate(lines)
                                            if re.match(r"^\[(\s*)system(\s*)]", iline)]
            idx = aux_find_end(lines, pos_system_start[ifile_name])
            for j in idx:
                pos_system_end[ifile_name].append(j)

            # [ molecules ]
            pos_molecules_start[ifile_name] = [index for index, iline in enumerate(lines)
                                               if re.match(r"^\[(\s*)molecules(\s*)]", iline)]
            idx = aux_find_end(lines, pos_molecules_start[ifile_name])
            for j in idx:
                pos_molecules_end[ifile_name].append(j)

            # [ moleculetype ]
            pos_moleculetype_start[ifile_name] = [index for index, iline in enumerate(lines)
                                                  if 'moleculetype' in iline.lower()]
            idx = aux_find_end(lines, pos_moleculetype_start[ifile_name])
            for j in idx:
                pos_moleculetype_end[ifile_name].append(j)

            # [ defaults ]
            pos_defaults_start[ifile_name] = [index for index, iline in enumerate(lines)
                                              if re.match(r"^\[(\s*)defaults(\s*)]", iline)]
            idx = aux_find_end(lines, pos_defaults_start[ifile_name])
            for j in idx:
                pos_defaults_end[ifile_name].append(j)

            # [ atomtypes]
            pos_atomtypes_start[ifile_name] = [index for index, iline in enumerate(lines)
                                               if re.match(r"^\[(\s*)atomtypes(\s*)]", iline)]
            idx = aux_find_end(lines, pos_atomtypes_start[ifile_name])
            for j in idx:
                pos_atomtypes_end[ifile_name].append(j)

            # [ bondtypes]
            pos_bondtypes_start[ifile_name] = [index for index, iline in enumerate(lines)
                                               if re.match(r"^\[(\s*)bondtypes(\s*)]", iline)]
            idx = aux_find_end(lines, pos_bondtypes_start[ifile_name])
            for j in idx:
                pos_bondtypes_end[ifile_name].append(j)

            # [pairtypes]
            pos_pairtypes_start[ifile_name] = [index for index, iline in enumerate(lines)
                                               if re.match(r"^\[(\s*)pairtypes(\s*)]", iline)]
            idx = aux_find_end(lines, pos_pairtypes_start[ifile_name])
            for j in idx:
                pos_pairtypes_end[ifile_name].append(j)

            # [angletypes]
            pos_angletypes_start[ifile_name] = [index for index, iline in enumerate(lines)
                                                if re.match(r"^\[(\s*)angletypes(\s*)]", iline)]
            idx = aux_find_end(lines, pos_angletypes_start[ifile_name])
            for j in idx:
                pos_angletypes_end[ifile_name].append(j)

            # [dihedraltypes]
            pos_dihedraltypes_start[ifile_name] = [index for index, iline in enumerate(lines)
                                                   if re.match(r"^\[(\s*)dihedraltypes(\s*)]", iline)]
            idx = aux_find_end(lines, pos_dihedraltypes_start[ifile_name])
            for j in idx:
                pos_dihedraltypes_end[ifile_name].append(j)

            # [ atoms ]
            ll = [index for index, iline in enumerate(lines) if re.match(r"^\[(\s*)atoms(\s*)]", iline)]
            if len(ll) > 0:
                pos_atoms_start[ifile_name] = ll
                idx = aux_find_end(lines, pos_atoms_start[ifile_name])
                for j in idx:
                    pos_atoms_end[ifile_name].append(j)

            # [ bonds ]
            ll = [index for index, iline in enumerate(lines) if re.match(r"^\[(\s*)bonds(\s*)]", iline)]
            if len(ll) > 0:
                pos_bonds_start[ifile_name] = ll
                idx = aux_find_end(lines, pos_bonds_start[ifile_name])
                for j in idx:
                    pos_bonds_end[ifile_name].append(j)

            # [ pairs ]
            ll = [index for index, iline in enumerate(lines) if re.match(r"^\[(\s*)pairs(\s*)]", iline)]
            if len(ll) > 0:
                pos_pairs_start[ifile_name] = ll
                idx = aux_find_end(lines, pos_pairs_start[ifile_name])
                for j in idx:
                    pos_pairs_end[ifile_name].append(j)

            # [ angles ]
            ll = [index for index, iline in enumerate(lines) if re.match(r"^\[(\s*)angles(\s*)]", iline)]
            if len(ll) > 0:
                pos_angles_start[ifile_name] = ll
                idx = aux_find_end(lines, pos_angles_start[ifile_name])
                for j in idx:
                    pos_angles_end[ifile_name].append(j)

            # [ dihedrals ]
            ll = [index for index, iline in enumerate(lines) if re.match(r"^\[(\s*)dihedrals(\s*)]", iline)]
            if len(ll) > 0:
                pos_dihedrals_start[ifile_name] = ll
                idx = aux_find_end(lines, pos_dihedrals_start[ifile_name])
                for j in idx:
                    pos_dihedrals_end[ifile_name].append(j)

    # [ defaults ]
    ntoken_defaults = aux_ntokens("[ defaults ]", pos_defaults_start)
    ntoken_atomtypes = aux_ntokens("[ atomtypes ]", pos_atomtypes_start)
    ntoken_bondtypes = aux_ntokens("[ bondtypes ]", pos_bondtypes_start)
    ntoken_angletypes = aux_ntokens("[ angletypes ]", pos_angletypes_start)
    ntoken_dihedraltypes = aux_ntokens("[ dihedraltypes ]", pos_dihedraltypes_start)
    ntoken_pairtypes = aux_ntokens("[ pairtypes ]", pos_pairtypes_start)
    ntoken_moleculetype = aux_ntokens("[ moleculetype ]", pos_moleculetype_start)
    if len(ntoken_bondtypes) != 0 or len(ntoken_angletypes) != 0 or\
       len(ntoken_dihedraltypes) != 0 or len(ntoken_pairtypes) != 0:
        isffbonded_terms = True
    else:
        isffbonded_terms = False

    # Generate forcefield.itp files
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m1 = "\n\t\tGenerating files in directory forcefield ({})\n".format(now)
    m = "\t\t" + len(m1) * "*"
    print(m1 + m) if logger is None else logger.info(m1 + m)

    with open(os.path.join(dir_target, "forcefield.itp"), 'w') as fff:

        try:
            idx_s = pos_defaults_start[ntoken_defaults[0]][0]
            idx_e = pos_defaults_end[ntoken_defaults[0]][0]
            with open(ntoken_defaults[0], "r") as finp:
                lines = finp.readlines()
                for jdx in range(idx_s, idx_e+1):
                    fff.writelines(lines[jdx])
        except IndexError:
            fff.writelines('[ defaults ]\n')
            fff.writelines('; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n')
            fff.writelines('1               3               no              1            0.5\n')

        fff.writelines("\n")
        fff.writelines('#include "ffnonbonded.itp"\n')
        if isffbonded_terms:
            fff.writelines('#include "ffbonded.itp"\n')

    # Generate ffnonbonded.itp files
    m1 = "\t\t\tGenerating ./forcefield/{}.itp".format("ffnonbonded")
    print(m1) if logger is None else logger.info(m1)
    with open(os.path.join(dir_target, "ffnonbonded.itp"), 'w') as fff:
        with open(ntoken_atomtypes[0], "r") as finp:
            lines = finp.readlines()
        label = lines[pos_atomtypes_start[ntoken_atomtypes[0]][0]]
        fff.writelines(label)
        for j in range(len(pos_atomtypes_start[ntoken_atomtypes[0]])):
            idx_s = pos_atomtypes_start[ntoken_atomtypes[0]][j]
            idx_e = pos_atomtypes_end[ntoken_atomtypes[0]][j]
            for jdx in range(idx_s+1, idx_e+1):
                iline = lines[jdx]
                try:
                    # [atomtypes]
                    # atom type; bonded type; atomic number; m (u); q (e); particle type; V
                    #  ; W  (bonded type and atomic number are optional)
                    # Sometimes the line containing atoms has comments at the end after a ";" symbol
                    idx_semicolon = iline.find(";")
                    if idx_semicolon != -1 and idx_semicolon == 0:
                        fff.writelines(iline)
                        continue
                    elif idx_semicolon != -1 and idx_semicolon != 0:
                        extraline = iline[idx_semicolon:]
                        iline = iline[:idx_semicolon]
                    else:
                        extraline = ""

                    nelements = len(iline.split())

                    if nelements == 6:
                        mass = float(iline.split()[1])
                        # element = [i for i in atomic_mass if mass + 0.1 > atomic_mass[i] > mass - 0.1][0]
                        e_itp = iline.split()[0]
                        if e_itp[0:2] != "CH":
                            element = [i for i in atomic_mass if np.isclose(mass, atomic_mass[i], atol=10 ** (-1))][0]
                        else:
                            element = "C"
                        at_number = atomic_number[element]
                        iline = '{0:12s} '.format(iline.split()[0]) + \
                                '{0:2s} '.format(element) + \
                                '{0:3d} '.format(int(at_number)) + \
                                '{0:8.4f} '.format(float(iline.split()[1])) + \
                                '{0:8.3f} '.format(float(iline.split()[2])) + \
                                '{0:2s} '.format(iline.split()[3]) + \
                                '{0:8.4f} '.format(float(iline.split()[4])) + \
                                '{0:8.4f}\n'.format(float(iline.split()[5]))

                    elif nelements == 7:
                        mass = float(iline.split()[2])
                        try:
                            int(iline.split()[1])
                            e_itp = iline.split()[0]
                            if e_itp[0:2] != "CH":
                                element = [i for i in atomic_mass if np.isclose(mass, atomic_mass[i], atol=10 ** (-1))][0]
                            else:
                                element = "C"
                            bondtypename = iline.split()[0]
                            iline = '{0:12s} '.format(iline.split()[0]) + \
                                    '{0:12s} '.format(bondtypename) + \
                                    '{0:3d} '.format(int(iline.split()[1])) + \
                                    '{0:8.4f} '.format(float(iline.split()[2])) + \
                                    '{0:8.3f} '.format(float(iline.split()[3])) + \
                                    '{0:2s} '.format(iline.split()[4]) + \
                                    '{0:8.4f} '.format(float(iline.split()[5])) + \
                                    '{0:8.4f}\n'.format(float(iline.split()[6]))
                        except ValueError:
                            # element = [i for i in atomic_mass if mass + 0.1 > atomic_mass[i] > mass - 0.1][0]
                            element = [i for i in atomic_mass if np.isclose(mass, atomic_mass[i], atol=10 ** (-1))][0]
                            e_itp = iline.split()[0]
                            if e_itp[0:2] != "CH":
                                element = [i for i in atomic_mass if np.isclose(mass, atomic_mass[i], atol=10 ** (-1))][0]
                            else:
                                element = "C"
                            at_number = atomic_number[element]

                            iline = '{0:12s} '.format(iline.split()[0]) + \
                                    '{0:2s} '.format(iline.split()[1]) + \
                                    '{0:3d} '.format(int(at_number)) + \
                                    '{0:8.4f} '.format(float(iline.split()[2])) + \
                                    '{0:8.3f} '.format(float(iline.split()[3])) + \
                                    '{0:2s} '.format(iline.split()[4]) + \
                                    '{0:8.4f} '.format(float(iline.split()[5])) + \
                                    '{0:8.4f}\n'.format(float(iline.split()[6]))

                    elif nelements == 8:
                        iline = '{0:12s} '.format(iline.split()[0]) + \
                                '{0:2s} '.format(iline.split()[1]) + \
                                '{0:3d} '.format(int(iline.split()[2])) + \
                                '{0:8.4f} '.format(float(iline.split()[3])) + \
                                '{0:8.3f} '.format(float(iline.split()[4])) + \
                                '{0:2s} '.format(iline.split()[5]) + \
                                '{0:8.4f} '.format(float(iline.split()[6])) + \
                                '{0:8.4f}\n'.format(float(iline.split()[7]))

                    else:
                        m = "\n\t\t ERROR. Atomtypes is not well-formed.\n"
                        m += "\t\t   {}\n".format(iline)
                        print(m) if logger is None else logger.error(m)
                        exit()

                except KeyError:
                    pass

                if len(extraline) != 0:
                    fff.writelines(iline[:-1]+"  "+extraline)
                else:
                    fff.writelines(iline)

    # Generate ffbonded.itp files
    if isffbonded_terms:
        m1 = "\t\t\tGenerating ./forcefield/{}.itp".format("ffbonded")
        print(m1) if logger is None else logger.info(m1)
        with open(os.path.join(dir_target, "ffbonded.itp"), 'w') as fff:
            # bondtypes
            with open(ntoken_bondtypes[0], "r") as finp:
                lines = finp.readlines()
            label = lines[pos_bondtypes_start[ntoken_bondtypes[0]][0]]
            fff.writelines(label)
            for j in range(len(pos_bondtypes_start[ntoken_bondtypes[0]])):
                idx_s = pos_bondtypes_start[ntoken_bondtypes[0]][j]
                idx_e = pos_bondtypes_end[ntoken_bondtypes[0]][j]
                for jdx in range(idx_s+1, idx_e+1):
                    fff.writelines(lines[jdx])
            fff.writelines("\n")

            # angletypes
            with open(ntoken_angletypes[0], "r") as finp:
                lines = finp.readlines()
            label = lines[pos_angletypes_start[ntoken_angletypes[0]][0]]
            fff.writelines(label)
            for j in range(len(pos_angletypes_start[ntoken_angletypes[0]])):
                idx_s = pos_angletypes_start[ntoken_angletypes[0]][j]
                idx_e = pos_angletypes_end[ntoken_angletypes[0]][j]
                for jdx in range(idx_s+1, idx_e+1):
                    fff.writelines(lines[jdx])
            fff.writelines("\n")

            # dihedraltypes
            with open(ntoken_dihedraltypes[0], "r") as finp:
                lines = finp.readlines()
            label = lines[pos_dihedraltypes_start[ntoken_dihedraltypes[0]][0]]
            fff.writelines(label)
            for j in range(len(pos_dihedraltypes_start[ntoken_dihedraltypes[0]])):
                idx_s = pos_dihedraltypes_start[ntoken_dihedraltypes[0]][j]
                idx_e = pos_dihedraltypes_end[ntoken_dihedraltypes[0]][j]
                for jdx in range(idx_s+1, idx_e+1):
                    fff.writelines(lines[jdx])
            fff.writelines("\n")

            # pairtypes
            with open(ntoken_pairtypes[0], "r") as finp:
                lines = finp.readlines()
            label = lines[pos_pairtypes_start[ntoken_pairtypes[0]][0]]
            fff.writelines(label)
            for j in range(len(pos_pairtypes_start[ntoken_pairtypes[0]])):
                idx_s = pos_pairtypes_start[ntoken_pairtypes[0]][j]
                idx_e = pos_pairtypes_end[ntoken_pairtypes[0]][j]
                for jdx in range(idx_s+1, idx_e+1):
                    fff.writelines(lines[jdx])
            fff.writelines("\n")
        itp_file = "./forcefield/ffbonded.itp"

    # Generate the itp file for each molecule type
    for imol in range(len(ntoken_moleculetype)):

        name_mol = os.path.splitext(os.path.split(ntoken_moleculetype[0])[-1])[0]
        m1 = "\t\t\tGenerating ./forcefield/{}.itp".format(name_mol)
        print(m1) if logger is None else logger.info(m1)

        with open(os.path.join(dir_target, name_mol + ".itp"), 'w') as fitp:
            line_header = "; Topology file write with create_library_monomer toolkit\n"
            line_header += ";                Javier Ramos (IEM-CSIC)\n\n"
            fitp.writelines(line_header)

            with open(ntoken_moleculetype[0], "r") as finp:
                lines = finp.readlines()
            # [ moleculetype ]
            label = lines[pos_moleculetype_start[ntoken_moleculetype[0]][0]]
            fitp.writelines(label)
            for j in range(len(pos_moleculetype_start[ntoken_moleculetype[0]])):
                idx_s = pos_moleculetype_start[ntoken_moleculetype[0]][j]
                idx_e = pos_moleculetype_end[ntoken_moleculetype[0]][j]
                for jdx in range(idx_s+1, idx_e+1):
                    fitp.writelines(lines[jdx])
            fitp.writelines("\n")

            # [ atoms ]
            label = lines[pos_atoms_start[ntoken_moleculetype[0]][0]]
            fitp.writelines(label)
            for j in range(len(pos_atoms_start[ntoken_moleculetype[0]])):
                idx_s = pos_atoms_start[ntoken_moleculetype[0]][j]
                idx_e = pos_atoms_end[ntoken_moleculetype[0]][j]
                for jdx in range(idx_s+1, idx_e+1):
                    fitp.writelines(lines[jdx])
            fitp.writelines("\n")

            # [ bonds ]
            label = lines[pos_bonds_start[ntoken_moleculetype[0]][0]]
            fitp.writelines(label)
            for j in range(len(pos_bonds_start[ntoken_moleculetype[0]])):
                idx_s = pos_bonds_start[ntoken_moleculetype[0]][j]
                idx_e = pos_bonds_end[ntoken_moleculetype[0]][j]
                for jdx in range(idx_s+1, idx_e+1):
                    fitp.writelines(lines[jdx])
            fitp.writelines("\n")

            # [ pairs ]
            try:
                label = lines[pos_pairs_start[ntoken_moleculetype[0]][0]]
                fitp.writelines(label)
                for j in range(len(pos_pairs_start[ntoken_moleculetype[0]])):
                    idx_s = pos_pairs_start[ntoken_moleculetype[0]][j]
                    idx_e = pos_pairs_end[ntoken_moleculetype[0]][j]
                    for jdx in range(idx_s+1, idx_e+1):
                        fitp.writelines(lines[jdx])
                fitp.writelines("\n")
                is_write_pairs = True
            except IndexError:  # No [pairs] section
                is_write_pairs = False
                pass

            # [ angles ]
            label = lines[pos_angles_start[ntoken_moleculetype[0]][0]]
            fitp.writelines(label)
            for j in range(len(pos_angles_start[ntoken_moleculetype[0]])):
                idx_s = pos_angles_start[ntoken_moleculetype[0]][j]
                idx_e = pos_angles_end[ntoken_moleculetype[0]][j]
                for jdx in range(idx_s+1, idx_e+1):
                    fitp.writelines(lines[jdx])
            fitp.writelines("\n")

            # [ dihedrals ]
            label = lines[pos_dihedrals_start[ntoken_moleculetype[0]][0]]
            fitp.writelines(label)
            for j in range(len(pos_dihedrals_start[ntoken_moleculetype[0]])):
                idx_s = pos_dihedrals_start[ntoken_moleculetype[0]][j]
                idx_e = pos_dihedrals_end[ntoken_moleculetype[0]][j]
                for jdx in range(idx_s+1, idx_e+1):
                    fitp.writelines(lines[jdx])
            fitp.writelines("\n")

        # Correct charges and residues in the generated itp
        # Remove repeated labels in the file
        itp_file = os.path.join(dir_target, name_mol + ".itp")
        correct_itp(itp_file, residuefile,
                    is_avg_charge=is_avg_charge, net_charge_res_zero=is_zero_charge_res, logger=logger)
        return is_write_pairs

# =================================================================================================
def expand_itpfile(itp_ffname):

    # Check for includes in itp_ffname file
    # and expand the includes in a itp or top temporal file.
    include_files = [itp_ffname]
    lines_out = ''
    nincludes = 0
    path_to_itp_ffname = os.path.split(itp_ffname)[0]
    with open(itp_ffname, 'r') as fitp:
        lines_top = fitp.readlines()
        for item in lines_top:
            if item.count("#include") > 0:
                nincludes += 1
                tmp_item = item.replace('"', '')
                ffname = os.path.join(path_to_itp_ffname, tmp_item.split()[-1])
                include_files.append(ffname)
    #             with open(ffname, 'r') as finclude:
    #                 l_tmp = finclude.readlines()
    #                 for j in l_tmp:
    #                     lines_out += j
    #         else:
    #             lines_out += item
    #
    # # Replace the generic bonds for bonds with parameters
    #
    #
    #
    #
    # tmp_topology_ffname = os.path.join(path_to_itp_ffname, "tmp_topology.itp")
    # with open(tmp_topology_ffname, 'w') as fout_tmp:
    #     fout_tmp.writelines(lines_out)

    return nincludes, include_files


# =================================================================================================
def check_pdb_itp_consistency(untyped_mol, pdb_filename, top_ffname, logger=None):
    """
    Check that the connectivity and the order of atoms is consistent
    between the pdb and itp files
    Returns:

    """

    # Check for includes in top_ffname file
    # and expand the includes in a itp or top temporal file.
    nincludes, include_files = expand_itpfile(top_ffname)
    itp_ffname = ""
    if nincludes > 1:
        for ifile in include_files:
            with open(ifile, 'r') as fitp:
                lines = fitp.readlines()
                pos_bonds = [index for index, iline in enumerate(lines) if re.match(r"^\[(\s*)bonds(\s*)]", iline)]
                if len(pos_bonds) != 0:
                    if itp_ffname != "":
                        m = "\n\t\t ERROR. Several itps contain bond information.\n"
                        m += "\t\t       ITP file 1: {}\n".format(itp_ffname)
                        m += "\t\t       ITP file 2: {}".format(ifile)
                        print(m) if logger is None else logger.error(m)
                        exit()
                    itp_ffname = ifile
    else:
        itp_ffname = include_files[0]

    # Check Connectivity
    pdb_graphdict = untyped_mol._topology._graphdict
    itp_graphdict = defaultdict(list)
    with open(itp_ffname, 'r') as fitp:
        lines = fitp.readlines()
        pos_bonds = [index for index, iline in enumerate(lines) if re.match(r"^\[(\s*)bonds(\s*)]", iline)]

        idx = pos_bonds[0] + 1
        while True:
            iline = lines[idx]
            if len(iline.replace(" ", "")) < 2 or lines[idx].find("[") != -1:
                break
            if re.match("^;", iline.strip()):
                idx += 1
                continue
            i, j, *_ = iline.split()
            itp_graphdict[int(i) - 1].append(int(j) - 1)
            itp_graphdict[int(j) - 1].append(int(i) - 1)
            idx += 1
        itp_graphdict = dict(sorted(itp_graphdict.items()))

        for key, item in itp_graphdict.items():
            itp_graphdict[key] = sorted(item)
        if not itp_graphdict == pdb_graphdict:
            m = "\n\t\t ERROR. PDB connection seems to be different to ITP connection.\n"
            m += "\t\t       ITP file: {}\n".format(itp_ffname)
            m += "\t\t       PDB file: {}".format(pdb_filename)
            print(m) if logger is None else logger.error(m)
            return False, include_files, itp_ffname

    # Check Connectivity
    pdb_elementsdict = untyped_mol._atom3d_element
    itp_elementsdict = defaultdict()
    with open(itp_ffname, 'r') as fitp:
        lines = fitp.readlines()
        pos_atoms = [index for index, iline in enumerate(lines) if re.match(r"^\[(\s*)atoms(\s*)]", iline)]

        idx = pos_atoms[0] + 1
        iat = 0
        while True:
            iline = lines[idx]
            if len(iline.replace(" ", "")) < 2 or lines[idx].find("[") != -1:
                break
            if re.match("^;", iline.strip()):
                idx += 1
                continue
            mass = float(iline.split()[7])
            e_itp = iline.split()[4]
            if e_itp[0] != "_":
                # element = [i for i in atomic_mass if atomic_mass[i] == mass][0]
                element = [i for i in atomic_mass if np.isclose(mass, atomic_mass[i], atol=10 ** (-1))][0]
            elif e_itp in ["_CH2", "_CH3", "_CH"]:
                element = "C"
            else:
                m = "\n\t\t ERROR. United atom {} does not exists.\n".format(e_itp)
                m += "\t\t       ITP file: {}\n".format(itp_ffname)
                m += "\t\t       PDB file: {}".format(pdb_filename)
                print(m) if logger is None else logger.error(m)
                # return False, include_files, itp_ffname
                exit()
            itp_elementsdict[iat] = element
            idx += 1
            iat += 1
    if not itp_elementsdict == pdb_elementsdict:
        m = "\n\t\t ERROR. PDB order for atoms seems to be  different to ITP file.\n"
        m += "\t\t       ITP file: {}\n".format(itp_ffname)
        m += "\t\t       PDB file: {}".format(pdb_filename)
        print(m) if logger is None else logger.error(m)
        # return False, include_files, itp_ffname
        exit()

    m = "\t\t\t PDB and ITP files seem to be consistent.\n"
    print(m) if logger is None else logger.info(m)

    return True, include_files, itp_ffname


# =================================================================================================
def correct_itp(itpfile_input, residue_file, is_avg_charge=True, net_charge_res_zero=False, logger=None):

    """
        # Correct charges and residues in the generated itp
        # Remove repeated labels in the file
    """

    def mean_and_std_dev(*lists):
        # Check if lists have the same size
        size = len(lists[0])
        if not all(len(lst) == size for lst in lists):
            raise ValueError("Lists must be of the same size")

        means = []
        std_devs = []

        # Calculate mean and standard deviation for each element
        for ii in range(size):
            elements = [lst[ii] for lst in lists]
            means.append(statistics.mean(elements))
            std_devs.append(statistics.stdev(elements))

        return means, std_devs

    # Get the residues for each atom
    residue_dict = defaultdict(list)
    if residue_file is not None:
        if os.path.splitext(residue_file)[-1] == ".pdb":
            with open(residue_file, "r") as fpdb:
                lines = fpdb.readlines()
                for iline in lines:
                    if iline.find("ATOM") == -1 and iline.find("HETATM") == -1:
                        continue
                    iat = int(iline.split()[1])
                    resname = iline.split()[3]
                    resnumber = int(iline.split()[5])
                    residue_dict[iat] = [resnumber, resname]
        else:
            m = "\n\t\t ERROR. Residue file must be PDB format.\n"
            print(m) if logger is None else logger.error(m)
            exit()

    # Compile residues and charges
    m = "\t\t\tAssigning new residues from PDB to ITP"
    print(m) if logger is None else logger.info(m)
    with open(itpfile_input, 'r') as fitp_inp:
        lines_itp = fitp_inp.readlines()

    charges_itp = defaultdict(list)
    residues_name = defaultdict(list)

    pos_atoms = [index for index, iline in enumerate(lines_itp)
                 if re.match(r"^\[(\s*)atoms(\s*)]", iline)][0]

    idx = pos_atoms + 1
    while True:
        iline = lines_itp[idx]
        if len(iline.replace(" ", "")) < 2 or lines[idx].find("[") != -1:
            break
        if re.match("^;", iline.strip()):
            idx += 1
            continue
        parts = iline.split()
        iat, itype, resnumber, resname, btype, cgr, charge, mass = parts[0:8]
        newresnumber = residue_dict[int(iat)][0]
        newresname = residue_dict[int(iat)][1]
        charges_itp[newresnumber].append(float(charge))
        if newresnumber not in residues_name[newresname]:
            residues_name[newresname].append(newresnumber)
        idx += 1

    m = "\t\t\tAveraging charges for residues of the same type: {}".format(is_avg_charge)
    print(m) if logger is None else logger.info(m)

    avg_charge_dict = defaultdict(list)
    avg_charge_dict_byatom = defaultdict()
    m = ""
    iat_idx = 0
    for iresname, values in residues_name.items():
        nelem = len(values)
        if nelem == 1:
            m += "\t\t\tThe charges for {} residues do not need to be averaged.\n" \
                .format(iresname)
            avg_charge_dict[iresname] = charges_itp[values[0]]
        else:
            m += "\t\t\tThe charges for {} residues are averaged. Number of residues: {}\n" \
                .format(iresname, nelem)
            ll = []
            for i in range(0, nelem):
                ll.append(charges_itp[values[i]])
            mean, std = mean_and_std_dev(*ll)
            avg_charge_dict[iresname] = mean
            for i, j in zip(mean, std):
                m += "\t\t\t\t{0:8.4f} +- {1:8.4f}\n".format(i, j)

    # If the net charge of the residues is not zero,
    # the excess of charge will be distribute among all atoms in the
    m1 = "\t\t\tRemove excess of charge in the residues of the {}: {}".\
        format(itpfile_input, net_charge_res_zero)
    print(m1) if logger is None else logger.info(m1)
    if net_charge_res_zero:
        for ires, charges in avg_charge_dict.items():
            original_sum = sum(charges)
            adjusment = original_sum / len(charges)
            adjustment_charges = [i - adjusment for i in charges]
            rounded_charges = [round(i, 4) for i in adjustment_charges]
            final_sum = round(sum(rounded_charges), 4)
            for i in range(len(rounded_charges)):
                rounded_charges[i] = round(rounded_charges[i] - final_sum / len(rounded_charges), 4)
            final_sum = round(sum(rounded_charges), 4)
            # Redistribute the excess of charge between the atoms
            if final_sum != 0.0:
                if final_sum < 0:
                    sign = 1
                else:
                    sign = -1
                iat = 0
                while round(np.abs(final_sum), 4) > 0.0:
                    rounded_charges[iat] = rounded_charges[iat] + sign*0.0001
                    final_sum = final_sum + sign*0.0001
                    iat += 1
            avg_charge_dict[ires] = rounded_charges

    for iresname, values in residues_name.items():
        nelem = len(values)
        for ielem in range(0, nelem):
            for icharge in avg_charge_dict[iresname]:
                avg_charge_dict_byatom[iat_idx] = icharge
                iat_idx += 1

    print(m) if logger is None else logger.info(m)

    # Change the residues and charges (write the file)
    with open(itpfile_input, 'w') as fitp_out:

        pos_atoms = [index for index, iline in enumerate(lines_itp)
                     if re.match(r"^\[(\s*)atoms(\s*)]", iline)][0]

        for i in range(0, pos_atoms + 1):
            fitp_out.writelines(lines_itp[i])

        idx = pos_atoms + 1
        qtotal = 0.0
        changeres = None
        while True:
            iline = lines_itp[idx]
            if len(iline.replace(" ", "")) < 2 or lines_itp[idx].find("[") != -1:
                break
            if re.match("^;", iline.strip()):
                fitp_out.writelines(lines_itp[idx])
                idx += 1
                continue
            parts = iline.split()
            iat, itype, resnumber, resname, btype, cgr, charge, mass = parts[0:8]
            newresnumber = residue_dict[int(iat)][0]
            newresname = residue_dict[int(iat)][1]
            if changeres is None or changeres != newresnumber:
                changeres = newresnumber
                fitp_out.writelines("; RES {}: {}\n".format(newresnumber, newresname))
            if is_avg_charge:
                newcharge = avg_charge_dict_byatom[int(iat)-1]
            else:
                newcharge = charge
            qtotal += newcharge
            e_itp = iline.split()[4]
            if e_itp[0] != "_":
                element = [i for i in atomic_mass if np.isclose(float(mass), atomic_mass[i], atol=10 ** (-1))][0]
            elif e_itp in ["_CH2", "_CH3", "_CH"]:
                element = "C"
            else:
                m = "\n\t\t ERROR. United atom {} does not exists.\n".format(e_itp)
                m += "\t\t       ITP file: {}\n".format(itpfile_input)
                print(m) if logger is None else logger.error(m)
                # return False, include_files, itp_ffname
                exit()

            newline = "{0:>6d}   {1:>10s}   {2:>4d}  {3:>5s}   {4:>5s}  " \
                      "{5:>3d}  {6:8.5f}  {7:8.4f} ;qt = {8:8.5f}\n".\
                format(int(iat), itype, newresnumber, newresname,
                       element, int(cgr), float(newcharge), float(mass), float(qtotal))
            fitp_out.writelines(newline)
            idx += 1

        for i in range(idx, len(lines_itp)):
            fitp_out.writelines(lines_itp[i])

    # The order of atoms in bonds section is very important.
    # The pair [i,j] must obey i<j
    with open(itpfile_input, 'r') as f2_inp:
        lines_itp = f2_inp.readlines()

    with open(itpfile_input, 'w') as fitp_out2:

        pos_bonds = [index for index, iline in enumerate(lines_itp)
                     if re.match(r"^\[(\s*)bonds(\s*)]", iline)][0]

        for i in range(0, pos_bonds + 1):
            fitp_out2.writelines(lines_itp[i])

        idx = pos_bonds + 1
        while True:
            iline = lines_itp[idx]
            if len(iline.replace(" ", "")) < 2 or lines_itp[idx].find("[") != -1:
                break
            if re.match("^;", iline.strip()):
                fitp_out2.writelines(lines_itp[idx])
                idx += 1
                continue
            iat, jat, *other = iline.split()
            iat = int(iat)
            jat = int(jat)
            if int(iat) > int(jat):
                newline = "{0:>8d} {1:>8d}   {2:s}\n".format(jat, iat, ' '.join(item for item in other))
            else:
                newline = "{0:>8d} {1:>8d}   {2:s}\n".format(iat, jat, ' '.join(item for item in other))
            fitp_out2.writelines(newline)
            idx += 1

        for i in range(idx, len(lines_itp)):
            fitp_out2.writelines(lines_itp[i])

    # Check and comments repeated keywords in the itp
    label_list = ["atoms", "moleculetype", "bonds", "angles", "dihedrals", "pairs"]
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m1 = "\n\t\tChecking for repeated labels ({})\n".format(now)
    m = "\t\t" + len(m1) * "*"
    print(m1 + m) if logger is None else logger.info(m1 + m)

    m = ""
    with open(itpfile_input, 'r') as fnewitp_out:
        lines_newitp = fnewitp_out.readlines()

        idx_not_write = []
        for ilabel in label_list:
            pos_labels = [index for index, iline in enumerate(lines_newitp)
                          if re.match(r"^\[(\s*){}(\s*)]".format(ilabel), iline)]
            if len(pos_labels) > 1:
                for j in range(1, len(pos_labels)):
                    m += "\t\t\tLabel  [{}] is repeated. Comment out all instances except the first.\n".format(ilabel)
                    idx_file = pos_labels[j]
                    lines_newitp[idx_file] = "; " + lines_newitp[idx_file]
                    for i in range(1, 3):
                        iline = lines_newitp[idx_file - i]
                        if len(iline.replace(" ", "")) < 2:
                            idx_not_write.append(idx_file-i)

    print(m) if logger is None else logger.info(m)

    # Delete blank lines between labels
    for item in idx_not_write:
        lines_newitp.pop(item)
    with open(itpfile_input, 'w') as fcorritp_out:
        fcorritp_out.writelines(lines_newitp)


# =================================================================================================
def correct_repeat_labels_itp(itpfile_input, logger=None):

    """
    Check and comment out repeated keywords in the ITP file.

    Parameters:
    itpfile_input (str): Path to the ITP file to be checked and corrected.
    label_list (list): List of labels to check for repetitions.
    logger (Logger, optional): Logger for logging messages. If None, print to console.

    Example: label_list = ["atoms", "moleculetype", "bonds", "angles", "dihedrals", "pairs"]
    """

    # Get current date and time
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

    # Prepare the initial log message
    m1 = "\n\t\tChecking for repeated labels ({})\n".format(now)
    m = "\t\t" + len(m1) * "*"
    print(m1 + m) if logger is None else logger.info(m1 + m)

    # Find labels in the file
    label_list = list()
    with open(itpfile_input, 'r') as fitp:
        for iline in fitp:
            iline.strip()
            if iline.startswith(";") and "[" in iline and "]" in iline:
                continue
            elif "[" in iline and "]" in iline:    # Check if the line contains the character '[' and ']'
                label_list.append(iline.strip())
            else:
                continue
        fitp.seek(0)
        lines_olditp = fitp.readlines()
        fitp.seek(0)
    label_list = sorted(list(set(label_list)))

    m = ""
    pos_labels_start = defaultdict(list)
    pos_labels_end = defaultdict(list)
    with open(itpfile_input, 'r') as fnewitp_out:
        lines_newitp = fnewitp_out.readlines()
        # Check for each label in the label list
        for ilabel in label_list:
            # Find positions of the current label in the file
            pos_labels_start[ilabel] = [index for index, iline in enumerate(lines_newitp)
                                        if ilabel == iline.strip()]
            for j in range(0, len(pos_labels_start[ilabel])):
                if len(pos_labels_start[ilabel]) > 1:
                    m += "\t\t\tLabel  [{}] is repeated. Comment out all instances except the first.\n".format(ilabel)
                idx = pos_labels_start[ilabel][j] + 1
                while True:
                    if idx >= len(lines_newitp):
                        break
                    iline = lines_newitp[idx]
                    if len(iline.replace(" ", "")) < 2:
                        break
                    idx += 1
                pos_labels_end[ilabel].append(idx-1)

    print(m) if logger is None else logger.info(m)

    # # Remove blank lines between labels
    # for item in idx_not_write:
    #     lines_newitp.pop(item)

    # Write the corrected content back to the ITP file
    lines_newitp = ""
    with open(itpfile_input, 'w') as fcorritp_out:
        for ilabel in label_list:
            lines_newitp += "{}\n".format(ilabel)
            for j in range(len(pos_labels_start[ilabel])):
                for idx_iline in range(pos_labels_start[ilabel][j]+1, pos_labels_end[ilabel][j] + 1):
                    lines_newitp += lines_olditp[idx_iline]
            lines_newitp += "\n"
        fcorritp_out.writelines(lines_newitp)
