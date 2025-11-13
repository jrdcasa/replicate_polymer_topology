import os
import re
import datetime
import numpy as np
from collections import namedtuple, defaultdict

# =============================================================================
def clean_lammps(gro_filename, top_filename,  ffname, create_lmp_mc_singlechain):
    """
    This function groups terms in the lmp and inp files.

    Returns:

    """

    lmp_filename_new = os.path.splitext(gro_filename)[0]+"_clean.lmp"
    inp_filename_new = os.path.splitext(gro_filename)[0]+"_clean.inp"

    # Detect if there are more than one moleculetype in the top file. If yes, then
    # the program create an unwrap top and gro file, called tmp_uwrp.top and tmp_uwrp.gro files
    with open(top_filename, 'r') as ftop:
        lines_top = ftop.readlines()
        nmoltypes = sum(1 for s in lines_top if re.match(r"^\[(\s*)moleculetype(\s*)\]", s))
        if nmoltypes > 1:
            _unwrap_top_gro_files(gro_filename, top_filename)
            top_filename = "tmp_join.top"


    # Open top_filename and look for pair styles
    d_indices_lines = defaultdict()
    with open(top_filename, 'r') as ftop:
        lines_top = ftop.readlines()
        ftop.seek(0)
        for num, line in enumerate(ftop, 1):
            if re.match(r"^\[(\s*)atomtypes(\s*)\]", line):
                d_indices_lines["atomtypes"] = num - 1
            if re.match(r"^\[(\s*)atoms(\s*)\]", line):
                d_indices_lines["atoms"] = num - 1
            if re.match(r"^\[(\s*)bonds(\s*)\]", line):
                d_indices_lines["bonds"] = num - 1
            if re.match(r"^\[(\s*)angles(\s*)\]", line):
                d_indices_lines["angles"] = num - 1
            if re.match(r"^\[(\s*)dihedrals(\s*)\]", line):
                d_indices_lines["dihedrals"] = num - 1
            if re.match(r"^\[(\s*)molecules(\s*)\]", line):
                d_indices_lines["molecules"] = num - 1
            if re.match(r"^\[(\s*)atomtypes(\s*)\]", line):
                d_indices_lines["molecules"] = num - 1
            if re.match(r"^\[(\s*)pairs(\s*)\]", line):
                d_indices_lines["pairs"] = num - 1

    # Get box dimensions from gro file
    with open(gro_filename, 'r') as fgro:
        lines_gro = fgro.readlines()
        fgro.seek(0)

    energy_terms, is_charged_system = \
        _write_lmp_data(lmp_filename_new, ffname, top_filename, lines_top, lines_gro, d_indices_lines)
    _write_lmp_inp(inp_filename_new, ffname, energy_terms,  is_charged_system)

    if create_lmp_mc_singlechain:
        lmp_filename_new_mc = os.path.splitext(gro_filename)[0] + "_clean.data"
        _write_lmp_data_mc(lmp_filename_new_mc, ffname, top_filename, lines_top, lines_gro, d_indices_lines)


# =============================================================================
def _write_lmp_data(lmp_filename_new, ffname, top_filename, lines_top, lines_gro, d_indices_lines):

    # Conversion factors
    kjmol_to_kcalmol = 0.239006
    nm_to_angstrom = 10

    natomtypes, atomtypes =\
        _parse_top_atomtypes_section(lines_top, d_indices_lines["atomtypes"], ffname)

    natoms_per_mol, atoms_info,  is_charged_system, atoms_kind_info =\
        _parse_top_atoms_section(lines_top, d_indices_lines["atoms"])

    nmols = _parse_top_molecules_section(lines_top, d_indices_lines["molecules"])
    nbonds_per_mol, bonds_type, bonds_info =\
        _parse_top_bonds_section(lines_top, d_indices_lines["bonds"], atoms_kind_info)
    nangles_per_mol, angles_type, angles_info =\
        _parse_top_angles_section(lines_top, d_indices_lines["angles"], atoms_kind_info)
    ndihedrals_per_mol, nimpropers_per_mol, dihedrals_type, impropers_type, dihedrals_info, impropers_info = \
        _parse_top_dihedrals_section(lines_top, d_indices_lines["dihedrals"], atoms_kind_info)
    # Check if there is new atom pairs under the section [pairs] of the top file
    # try:
    #     _parse_top_pairs_section(lines_top, d_indices_lines["pairs"], dihedrals_info)
    # except KeyError:
    #     pass

    coords, box_info_gro = _parse_gro_info(lines_gro)

    # Write lmp data
    with open(lmp_filename_new, "w") as flmp:
        flmp.writelines("LAMMPS data "+ffname+" "+top_filename+"\n")
        flmp.writelines("\n")
        flmp.writelines("{0:<8d}  atoms\n".format(natoms_per_mol*nmols))
        flmp.writelines("{0:<8d}  bonds\n".format(nbonds_per_mol*nmols))
        flmp.writelines("{0:<8d}  angles\n".format(nangles_per_mol * nmols))
        flmp.writelines("{0:<8d}  dihedrals\n".format(ndihedrals_per_mol * nmols))
        flmp.writelines("{0:<8d}  impropers\n".format(nimpropers_per_mol * nmols))
        flmp.writelines("\n")
        flmp.writelines("{0:<8d}  atom types\n".format(natomtypes))
        flmp.writelines("{0:<8d}  bond types\n".format(len(bonds_type)))
        flmp.writelines("{0:<8d}  angle types\n".format(len(angles_type)))
        flmp.writelines("{0:<8d}  dihedral types\n".format(len(dihedrals_type)))
        flmp.writelines("{0:<8d}  improper types\n".format(len(impropers_type)))
        flmp.writelines("\n")
        # From nanometers to anstroms
        flmp.writelines("{0:>11.5f} {1:>11.5f}   xlo xhi\n".format(box_info_gro[0][0]*10,
                                                                   box_info_gro[1][0]*10+box_info_gro[0][0]*10))
        flmp.writelines("{0:>11.5f} {1:>11.5f}   ylo yhi\n".format(box_info_gro[0][1]*10,
                                                                   box_info_gro[1][1]*10+box_info_gro[0][1]*10))
        flmp.writelines("{0:>11.5f} {1:>11.5f}   zlo zhi\n".format(box_info_gro[0][2]*10,
                                                                   box_info_gro[1][2]*10+box_info_gro[0][2]*10))
        flmp.writelines("{0:>11.5f} {1:>11.5f} {2:>11.5f} xy xz yz\n"
                        .format(box_info_gro[1][3]*10, box_info_gro[1][4]*10, box_info_gro[1][5]*10))
        # Masses
        flmp.writelines("\n")
        flmp.writelines("Masses\n")
        flmp.writelines("\n")
        for itype in atomtypes:
            flmp.writelines("{0:<8d} {1:>10.7}   # {2:s}\n".format(itype.id, itype.mass, itype.type))

        # Atoms
        flmp.writelines("\n")
        flmp.writelines("Atoms\n")
        flmp.writelines("\n")
        iat_global = 0
        for imol in range(nmols):
            for iat_local in atoms_info:
                idx = [at.id for at in atomtypes if at.type == iat_local.type][0]
                flmp.writelines("{0:<8d} {1:<6d} {2:<3d} {3:>8.4f} {4:10.4f} {5:10.4f} {6:10.4f} #{7:s} \n".
                                format(iat_global+1, imol+1, idx, iat_local.charge, coords[iat_global,0]*10,
                                       coords[iat_global,1]*10,coords[iat_global,2]*10, iat_local.type))
                iat_global += 1

        # Bonds
        if nbonds_per_mol != 0:
            flmp.writelines("\n")
            flmp.writelines("Bonds\n")
            flmp.writelines("\n")
            ib_global = 0
            for imol in range(nmols):
                for ib_local in bonds_info:
                    ib_global1 = ib_local[1] + imol * natoms_per_mol
                    ib_global2 = ib_local[2] + imol * natoms_per_mol
                    flmp.writelines("{0:<8d} {1:<6d} {2:<8d} {3:<8d} # {4:s} {5:s}\n".
                                    format(ib_global+1,ib_local[0],
                                           ib_global1,ib_global2, ib_local[-2], ib_local[-1]))
                    ib_global += 1

        # Angles
        if nangles_per_mol != 0:
            flmp.writelines("\n")
            flmp.writelines("Angles\n")
            flmp.writelines("\n")
            ia_global = 0
            for imol in range(nmols):
                for ia_local in angles_info:
                    ia_global1 = ia_local[1] + imol * natoms_per_mol
                    ia_global2 = ia_local[2] + imol * natoms_per_mol
                    ia_global3 = ia_local[3] + imol * natoms_per_mol
                    flmp.writelines("{0:<8d} {1:<6d} {2:<8d} {3:<8d} {4:<8d} # {5:s} {6:s} {7:s}\n".
                                    format(ia_global+1,ia_local[0],
                                           ia_global1, ia_global2, ia_global3,
                                           ia_local[-3], ia_local[-2], ia_local[-1]))
                    ia_global += 1

        # Dihedrals
        if ndihedrals_per_mol != 0:
            flmp.writelines("\n")
            flmp.writelines("Dihedrals\n")
            flmp.writelines("\n")
            id_global = 0
            for imol in range(nmols):
                for id_local in dihedrals_info:
                    id_global1 = id_local[1] + imol * natoms_per_mol
                    id_global2 = id_local[2] + imol * natoms_per_mol
                    id_global3 = id_local[3] + imol * natoms_per_mol
                    id_global4 = id_local[4] + imol * natoms_per_mol
                    flmp.writelines("{0:<8d} {1:<6d} {2:<8d} {3:<8d} {4:<8d} {5:<8d} # {6:s} {7:s} {8:s} {9:s}\n".
                                    format(id_global+1,id_local[0],
                                           id_global1, id_global2, id_global3, id_global4,
                                           id_local[-4], id_local[-3], id_local[-2], id_local[-1]))
                    id_global += 1

        # Impropers
        if nimpropers_per_mol != 0:
            flmp.writelines("\n")
            flmp.writelines("Impropers\n")
            flmp.writelines("\n")
            id_global = 0
            for imol in range(nmols):
                for id_local in impropers_info:
                    id_global1 = id_local[1] + imol * natoms_per_mol
                    id_global2 = id_local[2] + imol * natoms_per_mol
                    id_global3 = id_local[3] + imol * natoms_per_mol
                    id_global4 = id_local[4] + imol * natoms_per_mol
                    flmp.writelines("{0:<8d} {1:<6d} {2:<8d} {3:<8d} {4:<8d} {5:<8d} # {6:s} {7:s} {8:s} {9:s}\n".
                                    format(id_global+1,id_local[0],
                                           id_global1, id_global2, id_global3, id_global4,
                                           id_local[-4], id_local[-3], id_local[-2], id_local[-1]))
                    id_global += 1



        # Pair coefficients
        flmp.writelines("\n")
        flmp.writelines("Pair Coeffs\n")
        flmp.writelines("\n")
        for itype in atomtypes:
            flmp.writelines("{0:<4d} {1:s} {2:>8.4f} {3:>8.4f}   #{4:s}\n"
                            .format(itype.id, itype.nbpairs,
                                    itype.epsilon*kjmol_to_kcalmol,
                                    itype.sigma*nm_to_angstrom,
                                    itype.type))

        # Bond coefficients
        if nbonds_per_mol != 0:
            flmp.writelines("\n")
            flmp.writelines("Bond Coeffs\n")
            flmp.writelines("\n")
            for itype in bonds_type:
                if itype.type == "Harmonic":
                    flmp.writelines("{0:<4d} harmonic   {1:>8.4f} {2:>8.4f}   #{3:s} {4:s}\n"
                                    .format(itype.idx_bond,
                                            itype.params_list[1]*0.5*kjmol_to_kcalmol/(nm_to_angstrom*nm_to_angstrom),
                                            itype.params_list[0]*nm_to_angstrom,
                                            itype.type_A, itype.type_B ))

        # Angle coefficients
        if nangles_per_mol != 0:
            flmp.writelines("\n")
            flmp.writelines("Angle Coeffs\n")
            flmp.writelines("\n")
            for itype in angles_type:
                if itype.type == "Harmonic":
                    flmp.writelines("{0:<4d} harmonic   {1:>8.4f} {2:>8.4f}   #{3:s} {4:s} {5:s}\n"
                                    .format(itype.idx_angle,
                                            itype.params_list[0]*0.5*kjmol_to_kcalmol,
                                            itype.params_list[1],
                                            itype.type_A, itype.type_B, itype.type_C ))

        # Dihedral coefficients
        if ndihedrals_per_mol != 0:
            flmp.writelines("\n")
            flmp.writelines("Dihedral Coeffs\n")
            flmp.writelines("\n")
            for itype in dihedrals_type:
                if itype.type == "Ryckaert-Bellemans" or itype.type == "multi/harmonic":
                    flmp.writelines("{0:<4d} multi/harmonic   {1:>8.4f} {2:>8.4f} {3:>8.4f} {4:>8.4f} {6:>8.4f}  "
                                    "#{7:s} {8:s} {9:s} {10:s}\n"
                                    .format(itype.idx_dihedral,
                                            itype.params_list[0] * kjmol_to_kcalmol,
                                            -itype.params_list[1] * kjmol_to_kcalmol,
                                            itype.params_list[2] * kjmol_to_kcalmol,
                                            -itype.params_list[3] * kjmol_to_kcalmol,
                                            itype.params_list[4] * kjmol_to_kcalmol,
                                            0.000,
                                            itype.type_A, itype.type_B, itype.type_C, itype.type_D ))
                elif itype.type == "nharmonic":
                    params = dict()
                    nterms = len(itype.params_list)
                    k0 = itype.params_list[0]
                    params['k0'] = 2. * itype.params_list[0] + itype.params_list[1] + \
                                   itype.params_list[3] + 2. * itype.params_list[4] + \
                                   itype.params_list[5] + itype.params_list[7] + 2. * itype.params_list[8]
                    params['k1'] = -itype.params_list[1] + 3. * itype.params_list[3] \
                                   - 5. * itype.params_list[5] + 7. * itype.params_list[7]
                    params['k2'] = 2.* itype.params_list[2] - 8. * itype.params_list[4] \
                                   + 18. * itype.params_list[6] - 32. * itype.params_list[8]
                    params['k3'] = -4. * itype.params_list[3] + 20. * itype.params_list[5] - 56. * itype.params_list[7]
                    params['k4'] = 8. * itype.params_list[4] - 48. * itype.params_list[6] + 160. * itype.params_list[8]
                    params['k5'] = -16. * itype.params_list[5] + 112. * itype.params_list[7]
                    params['k6'] = 32. * itype.params_list[6] - 256. * itype.params_list[8]
                    params['k7'] = -64. * itype.params_list[7]
                    params['k8'] = 128. * itype.params_list[8]

                    line = "{0:<4d} nharmonic {1:2d}".format(itype.idx_dihedral, nterms)
                    for key, value in params.items():
                        line += " {0:>8.4f}".format(value/4.184)
                    itype.type_A, itype.type_B, itype.type_C, itype.type_D
                    line += "#{0:s} {1:s} {2:s} {3:s} \n".format(itype.type_A, itype.type_B, itype.type_C, itype.type_D)
                    flmp.writelines(line)

        # Improper coefficients
        if nimpropers_per_mol != 0:
            flmp.writelines("\n")
            flmp.writelines("Improper Coeffs\n")
            flmp.writelines("\n")
            for itype in impropers_type:
                if itype.type == "Harmonic":
                    flmp.writelines("{0:<4d} harmonic   {1:>8.4f} {2:>8.4f} "
                                    "#{3:s} {4:s} {5:s} {6:s}\n"
                                    .format(itype.idx_improper,
                                            itype.params_list[0] * 0.5 * kjmol_to_kcalmol,
                                            np.abs(itype.params_list[1]),
                                            itype.type_A, itype.type_B, itype.type_C, itype.type_D ))

        n_energy_terms = {"natoms": natoms_per_mol*nmols,
                          "nbonds": nbonds_per_mol*nmols,
                          "nangles": nangles_per_mol * nmols,
                          "ndihedrals": ndihedrals_per_mol * nmols,
                          "nimpropers": nimpropers_per_mol * nmols}

        return n_energy_terms, is_charged_system

# =============================================================================
def _write_lmp_data_mc(lmp_filename_new, ffname, top_filename, lines_top, lines_gro, d_indices_lines):

    # Conversion factors
    kjmol_to_kcalmol = 0.239006
    nm_to_angstrom = 10

    natomtypes, atomtypes =\
        _parse_top_atomtypes_section(lines_top, d_indices_lines["atomtypes"], ffname)

    natoms_per_mol, atoms_info,  is_charged_system, atoms_kind_info =\
        _parse_top_atoms_section(lines_top, d_indices_lines["atoms"])

    nmols = _parse_top_molecules_section(lines_top, d_indices_lines["molecules"])

    nbonds_per_mol, bonds_type, bonds_info =\
        _parse_top_bonds_section(lines_top, d_indices_lines["bonds"], atoms_kind_info)

    nangles_per_mol, angles_type, angles_info =\
        _parse_top_angles_section(lines_top, d_indices_lines["angles"], atoms_kind_info)

    ndihedrals_per_mol, nimpropers_per_mol, dihedrals_type, impropers_type, dihedrals_info, impropers_info = \
        _parse_top_dihedrals_section(lines_top, d_indices_lines["dihedrals"], atoms_kind_info)

    # Check if there is new atom pairs under the section [pairs] of the top file
    try:
        _parse_top_pairs_section(lines_top, d_indices_lines["pairs"], dihedrals_info)
    except KeyError:
        pass

    coords, box_info_gro = _parse_gro_info(lines_gro)

    # Write lmp data
    with open(lmp_filename_new, "w") as flmp:
        flmp.writelines("LAMMPS data "+ffname+" "+top_filename+"\n")
        flmp.writelines("\n")
        flmp.writelines("{0:<8d}  atoms\n".format(natoms_per_mol*nmols))
        flmp.writelines("{0:<8d}  bonds\n".format(nbonds_per_mol*nmols))
        flmp.writelines("{0:<8d}  angles\n".format(nangles_per_mol * nmols))
        flmp.writelines("{0:<8d}  dihedrals\n".format(ndihedrals_per_mol * nmols))
        flmp.writelines("{0:<8d}  impropers\n".format(nimpropers_per_mol * nmols))
        flmp.writelines("\n")
        flmp.writelines("{0:<8d}  atom types\n".format(natomtypes))
        flmp.writelines("{0:<8d}  bond types\n".format(len(bonds_type)))
        flmp.writelines("{0:<8d}  angle types\n".format(len(angles_type)))
        flmp.writelines("{0:<8d}  dihedral types\n".format(len(dihedrals_type)))
        flmp.writelines("{0:<8d}  improper types\n".format(len(impropers_type)))
        flmp.writelines("\n")
        # From nanometers to angstroms
        flmp.writelines("{0:>11.5f} {1:>11.5f}   xlo xhi\n".format(box_info_gro[0][0]*10,
                                                                   box_info_gro[1][0]*10+box_info_gro[0][0]*10))
        flmp.writelines("{0:>11.5f} {1:>11.5f}   ylo yhi\n".format(box_info_gro[0][1]*10,
                                                                   box_info_gro[1][1]*10+box_info_gro[0][1]*10))
        flmp.writelines("{0:>11.5f} {1:>11.5f}   zlo zhi\n".format(box_info_gro[0][2]*10,
                                                                   box_info_gro[1][2]*10+box_info_gro[0][2]*10))
        # Masses
        flmp.writelines("\n")
        flmp.writelines("Masses\n")
        flmp.writelines("\n")
        for itype in atomtypes:
            flmp.writelines("{0:<8d} {1:>10.7}   # {2:s}\n".format(itype.id, itype.mass, itype.type))

        # Atoms
        flmp.writelines("\n")
        flmp.writelines("Atoms\n")
        flmp.writelines("\n")
        iat_global = 0
        for imol in range(nmols):
            for iat_local in atoms_info:
                idx = [at.id for at in atomtypes if at.type == iat_local.type][0]
                flmp.writelines("{0:<8d} {1:<6d} {2:<3d} {3:>8.4f} {4:10.4f} {5:10.4f} {6:10.4f} #{7:s} \n".
                                format(iat_global+1, imol+1, idx, iat_local.charge, coords[iat_global,0]*10,
                                       coords[iat_global,1]*10,coords[iat_global,2]*10, iat_local.type))
                iat_global += 1

        # Bonds
        if nbonds_per_mol != 0:
            flmp.writelines("\n")
            flmp.writelines("Bonds\n")
            flmp.writelines("\n")
            ib_global = 0
            for imol in range(nmols):
                for ib_local in bonds_info:
                    ib_global1 = ib_local[1] + imol * natoms_per_mol
                    ib_global2 = ib_local[2] + imol * natoms_per_mol
                    flmp.writelines("{0:<8d} {1:<6d} {2:<8d} {3:<8d} # {4:s} {5:s}\n".
                                    format(ib_global+1,ib_local[0],
                                           ib_global1,ib_global2, ib_local[-2], ib_local[-1]))
                    ib_global += 1



        # Angles
        if nangles_per_mol != 0:
            flmp.writelines("\n")
            flmp.writelines("Angles\n")
            flmp.writelines("\n")
            ia_global = 0
            for imol in range(nmols):
                for ia_local in angles_info:
                    ia_global1 = ia_local[1] + imol * natoms_per_mol
                    ia_global2 = ia_local[2] + imol * natoms_per_mol
                    ia_global3 = ia_local[3] + imol * natoms_per_mol
                    flmp.writelines("{0:<8d} {1:<6d} {2:<8d} {3:<8d} {4:<8d} # {5:s} {6:s} {7:s}\n".
                                    format(ia_global+1,ia_local[0],
                                           ia_global1, ia_global2, ia_global3,
                                           ia_local[-3], ia_local[-2], ia_local[-1]))
                    ia_global += 1

        # Dihedrals
        if ndihedrals_per_mol != 0:
            flmp.writelines("\n")
            flmp.writelines("Dihedrals\n")
            flmp.writelines("\n")
            id_global = 0
            for imol in range(nmols):
                for id_local in dihedrals_info:
                    id_global1 = id_local[1] + imol * natoms_per_mol
                    id_global2 = id_local[2] + imol * natoms_per_mol
                    id_global3 = id_local[3] + imol * natoms_per_mol
                    id_global4 = id_local[4] + imol * natoms_per_mol
                    flmp.writelines("{0:<8d} {1:<6d} {2:<8d} {3:<8d} {4:<8d} {5:<8d} # {6:s} {7:s} {8:s} {9:s}\n".
                                    format(id_global+1,id_local[0],
                                           id_global1, id_global2, id_global3, id_global4,
                                           id_local[-4], id_local[-3], id_local[-2], id_local[-1]))
                    id_global += 1

        # Impropers
        if nimpropers_per_mol != 0:
            flmp.writelines("\n")
            flmp.writelines("Impropers\n")
            flmp.writelines("\n")
            id_global = 0
            for imol in range(nmols):
                for id_local in impropers_info:
                    id_global1 = id_local[1] + imol * natoms_per_mol
                    id_global2 = id_local[2] + imol * natoms_per_mol
                    id_global3 = id_local[3] + imol * natoms_per_mol
                    id_global4 = id_local[4] + imol * natoms_per_mol
                    flmp.writelines("{0:<8d} {1:<6d} {2:<8d} {3:<8d} {4:<8d} {5:<8d} # {6:s} {7:s} {8:s} {9:s}\n".
                                    format(id_global+1,id_local[0],
                                           id_global1, id_global2, id_global3, id_global4,
                                           id_local[-4], id_local[-3], id_local[-2], id_local[-1]))
                    id_global += 1



        ffname_b = os.path.splitext(os.path.split(ffname)[-1])[0].upper()

        # Generate all nb pairs
        geometric_nb_ff = ["OPLSAA", "LOPLSAA", "OPLSAA_AM"]
        arithmetic_nb_ff = ["TRAPPE-UA", "TRAPPE-UA_PETOXVAERD"]
        LJ_type_dict = {"OPLSAA": 1, "LOPLSAA": 1, "OPLSAA_AM":1, "TRAPPE-UA": 1, "TRAPPE-UA_PETOXVAERD": 1}
        try:
            LJ_type = LJ_type_dict[ffname_b]
        except KeyError:
            print("ERROR: prepare_lammps.py::_write_lmp_data_mc()")
            print("ERROR. {} force field has not LJ potential defined".format(ffname))
            print("ERROR. Revise LJ_type_dict")
            exit()

        lines_nb = []
        for i in range(0, natomtypes):
            for j in range(i, natomtypes):
                if i == j:
                    s = [s.sigma for s in atomtypes if s.id==i+1][0]
                    e = [e.epsilon*kjmol_to_kcalmol for e in atomtypes if e.id==i+1][0]
                else:
                    if ffname_b in geometric_nb_ff:
                        s = np.sqrt([s.sigma for s in atomtypes if s.id == i + 1][0]*
                                    [s.sigma for s in atomtypes if s.id == j + 1][0])
                        e = np.sqrt([e.epsilon*kjmol_to_kcalmol for e in atomtypes if e.id == i + 1][0]*
                                    [e.epsilon*kjmol_to_kcalmol for e in atomtypes if e.id == j + 1][0])

                    elif ffname_b in arithmetic_nb_ff:
                        s = 0.5*   ([s.sigma for s in atomtypes if s.id == i + 1][0]+
                                    [s.sigma for s in atomtypes if s.id == j + 1][0])
                        e = np.sqrt([e.epsilon*kjmol_to_kcalmol for e in atomtypes if e.id == i + 1][0]*
                                    [e.epsilon*kjmol_to_kcalmol for e in atomtypes if e.id == j + 1][0])
                    else:
                        print("ERROR: prepare_lammps.py::_write_lmp_data_mc()")
                        print("ERROR. {} force field does not mix rule defined".format(ffname))
                        print("ERROR. Revise geometric_nb_ff or arithmetic_nb_ff")
                        exit()
                t1 = [t1.type for t1 in atomtypes if t1.id == i+1][0]
                t2 = [t2.type for t2 in atomtypes if t2.id == j+1][0]
                lines_nb.append([i + 1, j + 1, e, s*10, LJ_type, t1, t2])

        # Pair coefficients
        flmp.writelines("\n")
        flmp.writelines("Pair Coeffs\n")
        flmp.writelines("\n")
        for iline_nb in lines_nb:

            flmp.writelines("{0:<4d} {1:<4d} {2:>8.4f} {3:>8.4f} {4:<2d} # {5:s} {6:s}\n"
                            .format(iline_nb[0], iline_nb[1], iline_nb[2], iline_nb[3],
                                    iline_nb[4], iline_nb[5], iline_nb[6]))

        # Bond coefficients
        if nbonds_per_mol != 0:
            flmp.writelines("\n")
            flmp.writelines("Bond Coeffs\n")
            flmp.writelines("\n")
            for itype in bonds_type:
                if itype.type == "Harmonic":
                    flmp.writelines("{0:<4d}  {1:>8.4f} {2:>8.4f}   #{3:s} {4:s}\n"
                                    .format(itype.idx_bond,
                                            itype.params_list[1]*0.5*kjmol_to_kcalmol/(nm_to_angstrom*nm_to_angstrom),
                                            itype.params_list[0]*nm_to_angstrom,
                                            itype.type_A, itype.type_B ))

        # Angle coefficients
        if nangles_per_mol != 0:
            flmp.writelines("\n")
            flmp.writelines("Angle Coeffs\n")
            flmp.writelines("\n")
            for itype in angles_type:
                if itype.type == "Harmonic":
                    flmp.writelines("{0:<4d}   {1:>8.4f} {2:>8.4f}   #{3:s} {4:s} {5:s}\n"
                                    .format(itype.idx_angle,
                                            itype.params_list[0]*0.5*kjmol_to_kcalmol,
                                            itype.params_list[1],
                                            itype.type_A, itype.type_B, itype.type_C ))

        # Dihedral coefficients
        if ndihedrals_per_mol != 0:
            flmp.writelines("\n")
            flmp.writelines("Dihedral Coeffs\n")
            flmp.writelines("\n")
            for itype in dihedrals_type:
                if itype.type == "Ryckaert-Bellemans" or itype.type == "multi/harmonic":
                    flmp.writelines("{0:<4d}  {1:>8.4f} {2:>8.4f} {3:>8.4f} {4:>8.4f} {6:>8.4f}  "
                                    "#{7:s} {8:s} {9:s} {10:s}\n"
                                    .format(itype.idx_dihedral,
                                            itype.params_list[0] * kjmol_to_kcalmol,
                                            -itype.params_list[1] * kjmol_to_kcalmol,
                                            itype.params_list[2] * kjmol_to_kcalmol,
                                            -itype.params_list[3] * kjmol_to_kcalmol,
                                            itype.params_list[4] * kjmol_to_kcalmol,
                                            0.000,
                                            itype.type_A, itype.type_B, itype.type_C, itype.type_D ))
                elif itype.type == "nharmonic":
                    params = dict()
                    nterms = len(itype.params_list)
                    k0 = itype.params_list[0]
                    params['k0'] = 2. * itype.params_list[0] + itype.params_list[1] + \
                                   itype.params_list[3] + 2. * itype.params_list[4] + \
                                   itype.params_list[5] + itype.params_list[7] + 2. * itype.params_list[8]
                    params['k1'] = -itype.params_list[1] + 3. * itype.params_list[3] \
                                   - 5. * itype.params_list[5] + 7. * itype.params_list[7]
                    params['k2'] = 2.* itype.params_list[2] - 8. * itype.params_list[4] \
                                   + 18. * itype.params_list[6] - 32. * itype.params_list[8]
                    params['k3'] = -4. * itype.params_list[3] + 20. * itype.params_list[5] - 56. * itype.params_list[7]
                    params['k4'] = 8. * itype.params_list[4] - 48. * itype.params_list[6] + 160. * itype.params_list[8]
                    params['k5'] = -16. * itype.params_list[5] + 112. * itype.params_list[7]
                    params['k6'] = 32. * itype.params_list[6] - 256. * itype.params_list[8]
                    params['k7'] = -64. * itype.params_list[7]
                    params['k8'] = 128. * itype.params_list[8]

                    line = "{0:<4d} nharmonic {1:2d}".format(itype.idx_dihedral, nterms)
                    for key, value in params.items():
                        line += " {0:>8.4f}".format(value/4.184)
                    line += "\n"
                    flmp.writelines(line)

        # Improper coefficients
        if nimpropers_per_mol != 0:
            flmp.writelines("\n")
            flmp.writelines("Improper Coeffs\n")
            flmp.writelines("\n")
            for itype in impropers_type:
                if itype.type == "Harmonic":
                    flmp.writelines("{0:<4d} harmonic   {1:>8.4f} {2:>8.4f} "
                                    "#{3:s} {4:s} {5:s} {6:s}\n"
                                    .format(itype.idx_improper,
                                            itype.params_list[0] * 0.5 * kjmol_to_kcalmol,
                                            np.abs(itype.params_list[1]),
                                            itype.type_A, itype.type_B, itype.type_C, itype.type_D ))

        n_energy_terms = {"natoms": natoms_per_mol*nmols,
                          "nbonds": nbonds_per_mol*nmols,
                          "nangles": nangles_per_mol * nmols,
                          "ndihedrals": ndihedrals_per_mol * nmols,
                          "nimpropers": nimpropers_per_mol * nmols}

        return n_energy_terms, is_charged_system


# =============================================================================
def _write_lmp_inp(inp_filename_new, ffname, energy_terms, ischargedsystem, units="real",
                   atom_style="full", temperature=450, timestep=1.0):

    with open(inp_filename_new, "w") as finp:
        finp.writelines("# LAMMPS Input Script created with replicate_polymer_topology\n")
        finp.writelines("# FF: {}\n".format(ffname))
        finp.writelines("\n")
        finp.writelines("#======1.Initialization======\n")
        finp.writelines("{0:<20s} {1:<40s}\n".format("units", units))
        finp.writelines("{0:<20s} {1:<40s}\n".format("boundary", "p p p"))
        finp.writelines("{0:<20s} {1:<40s}\n".format("dimension", "3"))
        finp.writelines("{0:<20s} {1:<40s}\n".format("atom_style", atom_style))
        finp.writelines("\n")

        # Get terms as a function of the force field
        ff = os.path.splitext(os.path.basename(ffname))[0]

        if ff.upper() == "OPLSAA" or ff.upper() == "OPLS" or ff.upper() == "LOPLSAA" or ff.upper() == "OPLSAA_AM" or ff.upper() == "LOPLSAA_AMESTER":
            ffterms = ["lj/cut/coul/long 10.0 10.0", "harmonic", "harmonic", "multi/harmonic", "harmonic"]
            ffterms_special_bonds = ["special_bonds lj 0.0 0.0 0.5 coul 0.0 0.0 0.5"]
            ffterms_pair_modify = ["pair_modify mix geometric #tail yes"]
            if ischargedsystem:
                ffterms_kspace = ["kspace_style pppm 1e-6\n #kspace_modify   diff  ad"]
            else:
                ffterms_kspace = ["kspace_style pppm 1e-6\nkspace_modify gewald 0.1"]
        elif ff.upper() == "TRAPPE-UA_PETOXVAERD" or ff.upper() == "TRAPPE-UA":
            ffterms = ["lj/cut/coul/long 10.0 10.0", "harmonic", "harmonic", "nharmonic multi/harmonic", "harmonic"]
            ffterms_special_bonds = ["special_bonds lj 0.0 0.0 0.0 coul 0.0 0.0 0.0"]
            ffterms_pair_modify = ["pair_modify mix arithmetic #tail yes"]
            if ischargedsystem:
                ffterms_kspace = ["kspace_style pppm 1e-6\n #kspace_modify   diff  ad"]
            else:
                ffterms_kspace = ["kspace_style pppm 1e-6\nkspace_modify gewald 0.1"]
        else:
            print("ERROR: prepare_lammps.py::_write_lmp_inpn()")
            print("ERROR. [bonds] section")
            print("ERROR. Force field {} is not yet implemented".format(ff.upper()))
            exit()

        finp.writelines("#======2.Force Field terms======\n")
        if energy_terms["natoms"] != 0:
            finp.writelines("pair_style     hybrid {}\n".format(ffterms[0]))
        else:
            finp.writelines("#pair_style     hybrid {}\n".format(ffterms[0]))
        if energy_terms["nbonds"] != 0:
            finp.writelines("bond_style     hybrid {}\n".format(ffterms[1]))
        else:
            finp.writelines("#bond_style     hybrid {}\n".format(ffterms[1]))
        if energy_terms["nangles"] != 0:
            finp.writelines("angle_style    hybrid {}\n".format(ffterms[2]))
        else:
            finp.writelines("#angle_style    hybrid {}\n".format(ffterms[2]))
        if energy_terms["ndihedrals"] != 0:
            finp.writelines("dihedral_style hybrid {}\n".format(ffterms[3]))
        else:
            finp.writelines("#dihedral_style hybrid {}\n".format(ffterms[3]))
        if energy_terms["nimpropers"] != 0:
            finp.writelines("improper_style hybrid {}\n".format(ffterms[4]))
        else:
            finp.writelines("#improper_style hybrid {}\n".format(ffterms[4]))
        finp.writelines("{}\n".format(ffterms_special_bonds[0]))
        finp.writelines("\n")

        finp.writelines("#======3.Data======\n")
        name = os.path.splitext(inp_filename_new)[0]+".lmp"
        finp.writelines("read_data {}\n".format(name))
        finp.writelines("\n")

        finp.writelines("#======4.Modify Pairs======\n")
        finp.writelines("{}\n".format(ffterms_pair_modify[0]))
        finp.writelines("{}\n".format(ffterms_kspace[0]))
        finp.writelines("\n")

        finp.writelines("#======5.Settings======\n")
        finp.writelines("neigh_modify every 10 delay 20 check yes\n")
        finp.writelines("thermo          {}\n".format(temperature))
        finp.writelines("thermo_style custom step cpu ebond eangle edihed eimp epair evdwl ecoul elong etail pe temp press lx ly lz density\n")
        finp.writelines("dump       1 all dcd      100  dump.dcd\n")
        finp.writelines("dump       2 all xyz    10000  dump_last.xyz\n")
        finp.writelines("dump       3 all atom   10000  dump_last.atom\n")
        finp.writelines("dump       4 all xtc      100  dump.xtc\n")
        finp.writelines("restart    10000 restart.data\n")
        finp.writelines("\n")

        finp.writelines("#======6.Minimization simulation======\n")
        finp.writelines("run 0\n")
        finp.writelines("#minimize 1.0e-8 1.0e-10 10000 10000\n")
        finp.writelines("\n")

        finp.writelines("#======6.NVT simulation======\n")
        finp.writelines("# timestep {}\n".format(timestep))
        finp.writelines("# fix 1 all nvt temp {} {} 100.0 drag 0.1\n".format(temperature, temperature))
        finp.writelines("# fix 2 all temp/csvr {} {} 500 54324  # Equivalent to v-rescale Gromacs\n".format(temperature, temperature))
        finp.writelines("# run 10000\n")
        finp.writelines("\n")

        finp.writelines("#======6.NPT simulation======\n")
        finp.writelines("# fix 1 all npt temp {} {} 100.0 iso 1.0 1.0 1000.0 drag 1.0\n".format(temperature, temperature))
        finp.writelines("# fix 2 all temp/csvr {} {} 500 54324  # Equivalent to v-rescale Gromacs\n".format(temperature, temperature))
        finp.writelines("# run 10000\n")
        finp.writelines("\n")

# =============================================================================
def _parse_top_atoms_section(lines_top, idx):

    natoms_per_mol = 0
    atoms_info_tuple = namedtuple("AtomInfo", ['idx_atom', 'type', 'resnumber', 'resname',
                                               'atomname', 'cgname', 'charge', 'mass'])
    atominfo = []
    charge_list = []
    atomkindinfo = []

    # Get lines with info ftom atoms
    idx += 1
    while True:
        iline = lines_top[idx]
        if len(iline.replace(" ", "")) < 2:
            break
        if not re.match("^;", iline.strip()) :
            tokens = iline.split()
            idx_atom = int(tokens[0])
            type = str(tokens[1])
            resnumber = int(tokens[2])
            resname = str(tokens[3])
            atomname = str(tokens[4])
            cgname = int(tokens[5])
            charge = float(tokens[6])
            charge_list.append(charge==0.0)
            mass = float(tokens[7])
            ainfo = atoms_info_tuple(idx_atom=idx_atom, type=type, resnumber=resnumber, resname=resname,
                                     atomname=atomname, cgname=cgname, charge=charge, mass=mass)
            atominfo.append(ainfo)
            atomkindinfo.append(ainfo.type)
            natoms_per_mol += 1
        idx += 1

    if any(charge_list):
        is_charged_system = False
    else:
        is_charged_system = True

    return natoms_per_mol, atominfo, is_charged_system, atomkindinfo

# =============================================================================
def _parse_top_bonds_section(lines_top, idx, atoms_kind_info):

    nbonds = 0
    bondtypes_tuple = namedtuple('Bond', ['idx_bond', 'type', 'params_list', 'type_A', 'type_B'])
    bondtypes = []
    bondinfo = []
    # Get lines with info ftom atoms
    idx += 1
    idxbond = 1
    while True:
        iline = lines_top[idx]
        if len(iline.replace(" ", "")) < 2:
            break
        if not re.match("^;", iline.strip()) :
            tokens = iline.split()
            # Harmonic bond potential
            if int(tokens[2]) == 1:
                idx1 = int(tokens[0])
                idx2 = int(tokens[1])
                k = float(tokens[3])
                b = float(tokens[4])
                # Types
                typeA = atoms_kind_info[idx1-1]
                typeB = atoms_kind_info[idx2-1]
                btype1 = bondtypes_tuple(idx_bond=idxbond, type="Harmonic", params_list=[k, b], type_A=typeA, type_B=typeB)
                btype2 = bondtypes_tuple(idx_bond=idxbond, type="Harmonic", params_list=[k, b], type_A=typeB, type_B=typeA)
                if not any(bb.type_A == btype1.type_A and bb.type_B == btype1.type_B for bb in bondtypes) and \
                   not any(bb.type_A == btype1.type_B and bb.type_B == btype1.type_A for bb in bondtypes):
                    bondtypes.append(btype1)
                    bondinfo.append([idxbond, idx1, idx2, typeA, typeB])
                    idxbond += 1
                else:
                    a1 = [bb.idx_bond for bb in bondtypes if bb.type_A==typeA and bb.type_B==typeB]
                    if len(a1) == 0:
                        a1 = [bb.idx_bond for bb in bondtypes if bb.type_A==typeB and bb.type_B==typeA]
                    bondinfo.append([a1[0], idx1, idx2, typeA, typeB])
            else:
                print("ERROR: prepare_lammps.py::_parse_top_bonds_section()")
                print("ERROR. [bonds] section")
                print("ERROR. Function {} is not yet implemented".format(int(tokens[2])))
                exit()
            nbonds += 1
        idx += 1

    return nbonds, bondtypes, bondinfo

# =============================================================================
def _parse_top_angles_section(lines_top, idx, atoms_kind_info):

    nangles = 0
    angletypes_tuple = namedtuple('Bend', ['idx_angle', 'type', 'params_list', 'type_A', 'type_B', 'type_C'])
    angletypes = []
    angleinfo = []

    # Get lines with info ftom atoms
    idx += 1
    idxangle = 1
    while True:
        iline = lines_top[idx]
        if len(iline.replace(" ", "")) < 2:
            break
        if not re.match("^;", iline.strip()) :
            tokens = iline.split()
            # Harmonic bend potential
            if int(tokens[3]) == 1:
                idx1 = int(tokens[0])
                idx2 = int(tokens[1])
                idx3 = int(tokens[2])
                ang = float(tokens[4])
                kbend = float(tokens[5])
                typeA = atoms_kind_info[idx1 - 1]
                typeB = atoms_kind_info[idx2 - 1]
                typeC = atoms_kind_info[idx3 - 1]
                atype1 = angletypes_tuple(idx_angle=idxangle, type="Harmonic", params_list=[kbend, ang],
                                         type_A=typeA, type_B=typeB, type_C=typeC)
                atype2 = angletypes_tuple(idx_angle=idxangle, type="Harmonic", params_list=[kbend, ang],
                                         type_A=typeC, type_B=typeB, type_C=typeA)

                if not any(aa.type_A == typeA and aa.type_B == typeB and aa.type_C == typeC for aa in angletypes) and \
                   not any(aa.type_A == typeC and aa.type_B == typeB and aa.type_C == typeA for aa in angletypes):
                    angletypes.append(atype1)
                    angleinfo.append([idxangle, idx1, idx2, idx3, typeA, typeB, typeC])
                    idxangle += 1
                else:
                    a1 = [aa.idx_angle for aa in angletypes if aa.type_A==typeA and
                                                               aa.type_B==typeB and
                                                               aa.type_C==typeC]
                    if len(a1) == 0:
                        a1 = [aa.idx_angle for aa in angletypes if aa.type_A == typeC and
                              aa.type_B == typeB and
                              aa.type_C == typeA]
                    angleinfo.append([a1[0], idx1, idx2, idx3, typeA, typeB, typeC])

            else:
                print("ERROR: prepare_lammps.py::_parse_top_angles_section()")
                print("ERROR. [angles] section")
                print("ERROR. Function {} is not yet implemented".format(int(tokens[3])))
                exit()

            nangles += 1
        idx += 1

    return nangles, angletypes, angleinfo

# =============================================================================
def _parse_top_dihedrals_section(lines_top, idx, atoms_kind_info):

    ndihedrals = 0
    nimpropers = 0
    dihedraltypes_tuple = namedtuple('Dihedral', ['idx_dihedral', 'type', 'params_list',
                                                  'type_A', 'type_B', 'type_C', 'type_D'])
    impropertypes_tuple = namedtuple('Improper', ['idx_improper', 'type', 'params_list',
                                                  'type_A', 'type_B', 'type_C', 'type_D'])
    dihedraltypes = []
    impropertypes = []
    dihedralinfo = []
    improperinfo = []

    # Get lines with info ftom atoms
    idx += 1
    idxdihedral = 1
    idximproper = 1
    while True:
        iline = lines_top[idx]
        if len(iline.replace(" ", "")) < 2:
            break
        if not re.match("^;", iline.strip()) :
            tokens = iline.split()
            if int(tokens[4]) == 3:           # Ryckaert-Bellemans dihedral for Trappe-UA or OPLS force fields.
                idx1 = int(tokens[0])
                idx2 = int(tokens[1])
                idx3 = int(tokens[2])
                idx4 = int(tokens[3])
                c0 = float(tokens[5])
                c1 = float(tokens[6])
                c2 = float(tokens[7])
                c3 = float(tokens[8])
                c4 = float(tokens[9])
                c5 = float(tokens[10])
                typeA = atoms_kind_info[idx1 - 1]
                typeB = atoms_kind_info[idx2 - 1]
                typeC = atoms_kind_info[idx3 - 1]
                typeD = atoms_kind_info[idx4 - 1]
                dtype1 = dihedraltypes_tuple(idx_dihedral=idxdihedral, type="Ryckaert-Bellemans",
                                             params_list=[c0, c1, c2, c3, c4, c5],
                                             type_A=typeA, type_B=typeB, type_C=typeC, type_D=typeD)
                dtype2 = dihedraltypes_tuple(idx_dihedral=idxdihedral, type="Ryckaert-Bellemans",
                                             params_list=[c0, c1, c2, c3, c4, c5],
                                             type_A=typeD, type_B=typeC, type_C=typeB, type_D=typeA)

                if not any(dd.type_A == typeA and dd.type_B == typeB and dd.type_C == typeC and dd.type_D == typeD
                           for dd in dihedraltypes) and \
                   not any(dd.type_A == typeD and dd.type_B == typeC and dd.type_C == typeB and dd.type_D == typeA
                                for dd in dihedraltypes):
                    dihedraltypes.append(dtype1)
                    dihedralinfo.append([idxdihedral, idx1, idx2, idx3, idx4, typeA, typeB, typeC, typeD])
                    idxdihedral += 1
                else:
                    d1 = [dd.idx_dihedral for dd in dihedraltypes if dd.type_A==typeA and
                                                                     dd.type_B==typeB and
                                                                     dd.type_C==typeC and
                                                                     dd.type_D==typeD]
                    if len(d1) == 0:
                        d1 = [dd.idx_dihedral for dd in dihedraltypes if dd.type_A == typeD and
                                                                         dd.type_B == typeC and
                                                                         dd.type_C == typeB and
                                                                         dd.type_D == typeA]
                    dihedralinfo.append([d1[0], idx1, idx2, idx3, idx4, typeA, typeB, typeC, typeD])

                ndihedrals += 1

            elif int(tokens[4]) == 9:   # Toxvaerd dihedral potential
                k = list()
                idx1 = int(tokens[0])
                idx2 = int(tokens[1])
                idx3 = int(tokens[2])
                idx4 = int(tokens[3])
                k.append(float(tokens[6]))
                typeA = atoms_kind_info[idx1 - 1]
                typeB = atoms_kind_info[idx2 - 1]
                typeC = atoms_kind_info[idx3 - 1]
                typeD = atoms_kind_info[idx4 - 1]
                label_previous = str(idx1) + "-" + str(idx2) + "-" + str(idx3) + "-" + str(idx4)
                jdx = idx + 1
                while True:
                    iline = lines_top[jdx]
                    tokens = iline.split()
                    if len(tokens) == 0:
                        idx = jdx - 1
                        break
                    try:
                        jdx1 = int(tokens[0])
                        jdx2 = int(tokens[1])
                        jdx3 = int(tokens[2])
                        jdx4 = int(tokens[3])
                    except ValueError:
                        break
                    label_current = str(jdx1) + "-" + str(jdx2) + "-" + str(jdx3) + "-" + str(jdx4)
                    if label_current != label_previous:
                        idx = jdx - 1
                        break
                    k.append(float(tokens[6]))
                    jdx += 1

                dtype1 = dihedraltypes_tuple(idx_dihedral=idxdihedral, type="nharmonic",
                                             params_list=k,
                                             type_A=typeA, type_B=typeB, type_C=typeC, type_D=typeD)
                dtype2 = dihedraltypes_tuple(idx_dihedral=idxdihedral, type="nharmonic",
                                             params_list=k,
                                             type_A=typeD, type_B=typeC, type_C=typeB, type_D=typeA)

                if not any(dd.type_A == typeA and dd.type_B == typeB and dd.type_C == typeC and dd.type_D == typeD
                           for dd in dihedraltypes) and \
                   not any(dd.type_A == typeD and dd.type_B == typeC and dd.type_C == typeB and dd.type_D == typeA
                                for dd in dihedraltypes):
                    dihedraltypes.append(dtype1)
                    dihedralinfo.append([idxdihedral, idx1, idx2, idx3, idx4, typeA, typeB, typeC, typeD])
                    idxdihedral += 1

                else:
                    d1 = [dd.idx_dihedral for dd in dihedraltypes if dd.type_A==typeA and
                                                                     dd.type_B==typeB and
                                                                     dd.type_C==typeC and
                                                                     dd.type_D==typeD]
                    if len(d1) == 0:
                        d1 = [dd.idx_dihedral for dd in dihedraltypes if dd.type_A == typeD and
                                                                         dd.type_B == typeC and
                                                                         dd.type_C == typeB and
                                                                         dd.type_D == typeA]
                    dihedralinfo.append([d1[0], idx1, idx2, idx3, idx4, typeA, typeB, typeC, typeD])

                ndihedrals += 1

            elif int(tokens[4]) == 2:                   # Improper angles
                k = list()
                idx1 = int(tokens[0])
                idx2 = int(tokens[1])
                idx3 = int(tokens[2])
                idx4 = int(tokens[3])
                a0 = float(tokens[5])
                k = float(tokens[6])
                typeA = atoms_kind_info[idx1 - 1]
                typeB = atoms_kind_info[idx2 - 1]
                typeC = atoms_kind_info[idx3 - 1]
                typeD = atoms_kind_info[idx4 - 1]
                dtype1 = impropertypes_tuple(idx_improper=idximproper, type="Harmonic",
                                             params_list=[k, a0],
                                             type_A=typeA, type_B=typeB, type_C=typeC, type_D=typeD)
                dtype2 = impropertypes_tuple(idx_improper=idximproper, type="Harmonic",
                                             params_list=[k, a0],
                                             type_A=typeD, type_B=typeC, type_C=typeB, type_D=typeA)

                if not any(dd.type_A == typeA and dd.type_B == typeB and dd.type_C == typeC and
                           dd.type_D == typeD and dd.params_list[1] == a0
                           for dd in impropertypes) and \
                   not any(dd.type_A == typeD and dd.type_B == typeC and dd.type_C == typeB and
                           dd.type_D == typeA and dd.params_list[1] == a0
                                for dd in impropertypes):
                    impropertypes.append(dtype1)
                    improperinfo.append([idximproper, idx1, idx2, idx3, idx4, typeA, typeB, typeC, typeD])
                    idximproper += 1
                else:
                    d1 = [dd.idx_improper for dd in impropertypes if dd.type_A==typeA and
                                                                     dd.type_B==typeB and
                                                                     dd.type_C==typeC and
                                                                     dd.type_D==typeD]
                    if len(d1) == 0:
                        d1 = [dd.idx_improper for dd in dihedraltypes if dd.type_A == typeD and
                                                                         dd.type_B == typeC and
                                                                         dd.type_C == typeB and
                                                                         dd.type_D == typeA]
                    improperinfo.append([d1[0], idx1, idx2, idx3, idx4, typeA, typeB, typeC, typeD])


                nimpropers += 1

            else:
                print("ERROR: prepare_lammps.py::_parse_top_dihedrals_section()")
                print("ERROR. [dihedrals] section")
                print("ERROR. Function {} is not yet implemented".format(int(tokens[4])))
                exit()

        idx += 1

    # for i in dihedraltypes:
    #     print(i)

    return ndihedrals, nimpropers, dihedraltypes, impropertypes, dihedralinfo, improperinfo

# =============================================================================
def _parse_gro_info(lines_gro):

    title = lines_gro[0]
    natoms = int(lines_gro[1])
    coords = np.zeros([natoms, 3])

    # Get coordinates
    for iatom in range(natoms):
        iline = lines_gro[iatom+2]
        coords[iatom, 0] = float(iline[20:28])
        coords[iatom, 1] = float(iline[28:36])
        coords[iatom, 2] = float(iline[36:44])

    # Max and minimum coordinates in all directions
    minx = np.min(coords[:, 0])
    miny = np.min(coords[:, 1])
    minz = np.min(coords[:, 2])
    maxx = np.max(coords[:, 0])
    maxy = np.max(coords[:, 1])
    maxz = np.max(coords[:, 2])

    # xx, yy, zz,
    box_line = lines_gro[-1].split()
    xx = float(box_line[0])
    yy = float(box_line[1])
    zz = float(box_line[2])
    try:
        xy = float(box_line[5])
        xz = float(box_line[7])
        yz = float(box_line[8])
    except:
        xy = 0.0
        xz = 0.0
        yz = 0.0

    # minx, miny, minz, maxx, maxy, maxz AND xx yy zz .. .. ..
    box_info_gro = [[minx, miny, minz, maxx, maxy, maxz], [xx, yy, zz, xy, xz, yz] ]

    return coords, box_info_gro

# =============================================================================
def _parse_top_molecules_section(lines_top, idx):

    nmols = 0

    # Get lines with info ftom atoms
    idx += 1
    while True:
        try:
            iline = lines_top[idx]
        except IndexError:
            break
        if len(iline.replace(" ", "")) < 2:
            break
        if not re.match("^;", iline.strip()):
            tokens = iline.split()
            nmols += int(tokens[1])
        idx += 1

    return nmols

# =============================================================================
def _parse_top_atomtypes_section(lines_top, idx, ffname):

    natomtypes = 0
    atomtypes_tuple = namedtuple("AtomTypes",['id', 'type', 'at_num', 'nbpairs', 'mass', 'charge', 'sigma', 'epsilon'])
    atomtypes = []

    # Get terms as a function of the force field
    ff = os.path.splitext(os.path.basename(ffname))[0]

    if ff.upper() == "OPLSAA" or ff.upper() == "OPLS" or ff.upper() == "LOPLSAA" or\
       ff.upper() == "TRAPPE-UA_PETOXVAERD" or ff.upper() == "TRAPPE-UA" or ff.upper() == "OPLSAA_AM"  or ff.upper() == "LOPLSAA_AMESTER":
        nbpairs = "lj/cut/coul/long"
    else:
        print("ERROR: prepare_lammps.py::_parse_top_atomtypes_section()")
        print("ERROR. [atomtypes] section")
        print("ERROR. Force field {} is not yet implemented".format(ff.upper()))
        exit()

    # Get lines with info from atoms
    idx += 1
    while True:
        iline = lines_top[idx]
        if len(iline.replace(" ", "")) < 2:
            break
        if not re.match("^;", iline.strip()) :
            tokens = iline.split()
            if len(tokens) != 7:
                print("ERROR: prepare_lammps.py::_parse_top_atomtypes_section()")
                print("ERROR. [atomtypes] section must have 7 elements. {} elements found".format(len(tokens)))
                exit()
            else:
                b = atomtypes_tuple(natomtypes+1, str(tokens[0]), int(tokens[1]), str(nbpairs),
                                    float(tokens[2]), float(tokens[3]),
                                    float(tokens[5]), float(tokens[6]))
                atomtypes.append(b)
            natomtypes += 1
        idx += 1

    return natomtypes, atomtypes

# =============================================================================
def _parse_top_pairs_section(lines_top, idx, dihedrals_info):

    idx += 1
    while True:
        iline = lines_top[idx]
        if len(iline.replace(" ", "")) < 2:
            break
        if not re.match("^;", iline.strip()) :
            tokens = iline.split()
            at1 = int(tokens[0])
            at4 = int(tokens[1])
            if [i[1] == at1 and i[4] == at4 for i in dihedrals_info] or \
               [i[1] == at4 and i[4] == at1 for i in dihedrals_info] :
                pass
            else:
                print("WARNING: prepare_lammps.py::_parse_top_pairs_section()")
                print("WARNING: A non-bonded pair (not 1-4) is defined in the [pairs] section, \n"
                      "         which is not included in the lammps data file. ")
                print("WARNING: Check and correct the problem.")
                print("WARNING: Introduce manually the information.")
                print("WARNING: {}".format(iline))
        idx += 1


# =============================================================================
def _unwrap_top_gro_files(gro_filename, top_filename):

    """
    This subroutine join the molecules in a top file in only one to make easy the transition from
    GROMACS to LAMMPS

    Args:
        gro_filename:
        top_filename:

    Returns:

    """

    uwrp_indices_lines = defaultdict(list)
    dict_moltypes = defaultdict()
    dict_typesmol = defaultdict()
    dict_atoms = defaultdict(list)
    dict_bonds = defaultdict(list)
    dict_pairs = defaultdict(list)
    dict_angles = defaultdict(list)
    dict_dihedrals = defaultdict(list)
    new_top_lines = []
    with open(top_filename, 'r') as ftop:
        lines_oldtop = ftop.readlines()
        ftop.seek(0)
        for num, line in enumerate(ftop, 1):
            if re.match(r"^\[(\s*)defaults(\s*)\]", line):
                uwrp_indices_lines["defaults"].append(num - 1)
            if re.match(r"^\[(\s*)atomtypes(\s*)\]", line):
                uwrp_indices_lines["atomtypes"].append(num - 1)
            if re.match(r"^\[(\s*)moleculetype(\s*)\]", line):
                uwrp_indices_lines["moleculetype"].append(num - 1)
            if re.match(r"^\[(\s*)atoms(\s*)\]", line):
                uwrp_indices_lines["atoms"].append(num - 1)
            if re.match(r"^\[(\s*)bonds(\s*)\]", line):
                uwrp_indices_lines["bonds"].append(num - 1)
            if re.match(r"^\[(\s*)angles(\s*)\]", line):
                uwrp_indices_lines["angles"].append(num - 1)
            if re.match(r"^\[(\s*)dihedrals(\s*)\]", line):
                uwrp_indices_lines["dihedrals"].append(num - 1)
            if re.match(r"^\[(\s*)pairs(\s*)\]", line):
                uwrp_indices_lines["pairs"].append(num - 1)
            if re.match(r"^\[(\s*)system(\s*)\]", line):
                uwrp_indices_lines["system"].append(num - 1)
            if re.match(r"^\[(\s*)molecules(\s*)\]", line):
                uwrp_indices_lines["molecules"].append(num - 1)

        # Read moleculetype
        imol = 0
        for imol_idx in uwrp_indices_lines["moleculetype"]:
            i = imol_idx+1
            while True:
                if len(lines_oldtop[i].strip()) < 2:
                    break
                if re.match(r"^;", lines_oldtop[i]):
                    i += 1
                    continue
                #dict_moltypes[lines_oldtop[i].split()[0]] = imol
                dict_typesmol[imol] = lines_oldtop[i].split()[0]
                imol += 1
                i += 1
        # Read atoms
        imol = 0
        for idx in uwrp_indices_lines["atoms"]:
            i = idx + 1
            while True:
                if len(lines_oldtop[i].strip()) < 2:
                    break
                if re.match(r"^;", lines_oldtop[i]):
                    i += 1
                    continue
                dict_atoms[dict_typesmol[imol]].append(lines_oldtop[i])
                i += 1
            imol += 1
        # Read bonds
        imol = 0
        for idx in uwrp_indices_lines["bonds"]:
            i = idx + 1
            while True:
                if len(lines_oldtop[i].strip()) < 2:
                    break
                if re.match(r"^;", lines_oldtop[i]):
                    i += 1
                    continue
                dict_bonds[dict_typesmol[imol]].append(lines_oldtop[i])
                i += 1
            imol += 1
        # Read pairs
        imol = 0
        for idx in uwrp_indices_lines["pairs"]:
            i = idx + 1
            while True:
                if len(lines_oldtop[i].strip()) < 2:
                    break
                if re.match(r"^;", lines_oldtop[i]):
                    i += 1
                    continue
                dict_pairs[dict_typesmol[imol]].append(lines_oldtop[i])
                i += 1
            imol += 1
        # Read angles
        imol = 0
        for idx in uwrp_indices_lines["angles"]:
            i = idx + 1
            while True:
                if len(lines_oldtop[i].strip()) < 2:
                    break
                if re.match(r"^;", lines_oldtop[i]):
                    i += 1
                    continue
                dict_angles[dict_typesmol[imol]].append(lines_oldtop[i])
                i += 1
            imol += 1

        # Read dihedrals
        imol = 0
        for idx in uwrp_indices_lines["dihedrals"]:
            i = idx + 1
            while True:
                if len(lines_oldtop[i].strip()) < 2:
                    break
                if re.match(r"^;", lines_oldtop[i].strip()):
                    i += 1
                    continue
                dict_dihedrals[dict_typesmol[imol]].append(lines_oldtop[i])
                i += 1
            imol += 1

        # Read [ molecules ] information
        imol = 0
        for idx in uwrp_indices_lines["molecules"]:
            i = idx + 1
            while True:
                try:
                    if len(lines_oldtop[i].strip()) < 2:
                        break
                except IndexError:
                    break
                if re.match(r"^;", lines_oldtop[i].strip()):
                    i += 1
                    continue
                namemol, numbermol = lines_oldtop[i].split()[0:2]
                numbermol = int(numbermol)
                nat_mol = len(dict_atoms[namemol])
                dict_moltypes[imol] = [namemol, numbermol, nat_mol]
                i += 1
                imol += 1
        # Precalculate map global atom index to molecule
        map_globalatidx_to_molname = defaultdict()
        iat_global = 0
        for key, values in dict_moltypes.items():
            namemol = values[0]
            numbermols = values[1]
            nat_mol = len(dict_atoms[namemol])
            for imol in range(0, numbermols):
                for _ in range(0, nat_mol):
                    map_globalatidx_to_molname[iat_global] = namemol
                    iat_global += 1

    # Get unwrapped data
    idx_atom_local = 0
    p_namemol = ""
    list_new_atoms = []
    local_map_list = []
    map_local_global = []
    map_order_molecules = []
    with open(gro_filename, 'r') as fgro:
        title = fgro.readline()
        natoms = int(fgro.readline())
        imol = 0
        idx_atom_global = 0
        for imoltype, values in dict_moltypes.items():
            namemol = values[0]
            nmols_imoltype = values[1]
            nat_imoltype = values[2]
            for inmol in range(0, nmols_imoltype):
                local_map_list = []
                for inat in range(0, nat_imoltype):
                    local_map_list.append(idx_atom_global)
                    iline = dict_atoms[namemol][inat]
                    tokens_iline = iline.split()
                    new_line = "{0:6d} ".format(idx_atom_global+1)
                    tokens_iline[2] = str(imol)
                    tokens_iline[3] = 'mol'
                    tokens_iline[5] = str(imol)
                    for itoken in tokens_iline[1:]:
                        new_line += itoken + " "
                    new_line += "\n"
                    list_new_atoms.append(new_line)
                    idx_atom_global += 1
                map_local_global.append(local_map_list)
                map_order_molecules.append(namemol)

    # print("dict_bonds", dict_bonds)
    # print(len(dict_bonds))
    # print("map_order_molecules", map_order_molecules)
    # print(len(map_order_molecules))
    # print("map_local_global", map_local_global)
    # print(len(map_local_global))
    # print("dict_moltypes", dict_moltypes)
    # print(len(dict_moltypes))

    # Write new top file
    with open("tmp_join.top", 'w') as fnewtop:
        # [defaults]
        idx = uwrp_indices_lines["defaults"][0]
        while True:
            iline = lines_oldtop[idx]
            if len(iline.strip()) < 2:
                break
            fnewtop.writelines(iline)
            idx += 1
        fnewtop.writelines("\n")
        # [atomtypes]
        idx = uwrp_indices_lines["atomtypes"][0]
        while True:
            iline = lines_oldtop[idx]
            if len(iline.strip()) < 2:
                break
            fnewtop.writelines(iline)
            idx += 1
        fnewtop.writelines("\n")
        # [moleculetype]
        idx = uwrp_indices_lines["moleculetype"][0]
        ll = []
        while True:
            iline = lines_oldtop[idx]
            if len(iline.strip()) < 2:
                break
            ll.append(iline)
            idx += 1
        for j in ll[:-1]:
            fnewtop.writelines(j)
        fnewtop.writelines("mol    {}\n".format(ll[-1].split()[-1]))
        fnewtop.writelines("\n")

        # [atoms]
        fnewtop.writelines("[ atoms ]\n")
        for iline in list_new_atoms:
            fnewtop.writelines(iline)
        fnewtop.writelines("\n")

        # [bonds]
        fnewtop.writelines("[ bonds ]\n")
        imol_count = 0
        for imol in map_order_molecules:
            for ibond in dict_bonds[imol]:
                tokens = ibond.split()
                iat1_local = int(tokens[0])
                iat2_local = int(tokens[1])
                iline = "{0:5d} {1:5d}   ".format(map_local_global[imol_count][iat1_local-1]+1,
                                                  map_local_global[imol_count][iat2_local-1]+1)
                for itoken in tokens[2:]:
                    iline += itoken+"   "
                iline += "\n"
                fnewtop.writelines(iline)
            imol_count += 1
        fnewtop.writelines("\n")

        # [pairs]
        fnewtop.writelines("[ pairs ]\n")
        imol_count = 0
        for imol in map_order_molecules:
            for ipair in dict_pairs[imol]:
                tokens = ipair.split()
                iat1_local = int(tokens[0])
                iat2_local = int(tokens[1])
                iline = "{0:5d} {1:5d}   ".format(map_local_global[imol_count][iat1_local-1]+1,
                                                  map_local_global[imol_count][iat2_local-1]+1)
                for itoken in tokens[2:]:
                    iline += itoken+"   "
                iline += "\n"
                fnewtop.writelines(iline)
            imol_count += 1
        fnewtop.writelines("\n")

        # [angles]
        fnewtop.writelines("[ angles ]\n")
        imol_count = 0
        for imol in map_order_molecules:
            for iangle in dict_angles[imol]:
                tokens = iangle.split()
                iat1_local = int(tokens[0])
                iat2_local = int(tokens[1])
                iat3_local = int(tokens[2])
                iline = "{0:5d} {1:5d}  {2:5d}   ".format(map_local_global[imol_count][iat1_local-1]+1,
                                                          map_local_global[imol_count][iat2_local-1]+1,
                                                          map_local_global[imol_count][iat3_local-1]+1)
                for itoken in tokens[3:]:
                    iline += itoken+"   "
                iline += "\n"
                fnewtop.writelines(iline)
            imol_count += 1
        fnewtop.writelines("\n")

        # [dihedrals]
        fnewtop.writelines("[ dihedrals ]\n")
        imol_count = 0
        for imol in map_order_molecules:
            for idih in dict_dihedrals[imol]:
                tokens = idih.split()
                iat1_local = int(tokens[0])
                iat2_local = int(tokens[1])
                iat3_local = int(tokens[2])
                iat4_local = int(tokens[3])
                iline = "{0:5d}  {1:5d}  {2:5d}   {3:5d}   ".format(map_local_global[imol_count][iat1_local-1]+1,
                                                                    map_local_global[imol_count][iat2_local-1]+1,
                                                                    map_local_global[imol_count][iat3_local-1]+1,
                                                                    map_local_global[imol_count][iat4_local-1]+1)
                for itoken in tokens[4:]:
                    iline += itoken+"   "
                iline += "\n"
                fnewtop.writelines(iline)
            imol_count += 1
        fnewtop.writelines("\n")

        fnewtop.writelines("[ system ]\n")
        fnewtop.writelines("Generic name\n")
        fnewtop.writelines("\n")

        fnewtop.writelines("[ molecules ]\n")
        fnewtop.writelines("mol 1\n")
        fnewtop.writelines("\n")


