import os
import re
import datetime
import numpy as np
from collections import namedtuple, defaultdict
from intermol.convert import _load_gromacs
from intermol.lammps.lammps_parser import LammpsParser


# =============================================================================
def prepare_lammps(basefile, title=None, nonbonded_style_defaults=None,
                   logger=None):

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m1 = "\t\tGenerating lmp and data files for {} ({})\n".format("LAMMPS", now)
    m = "\t\t" + len(m1) * "*" + "\n"
    print(m1 + m) if logger is None else logger.info(m1 + m)

    if nonbonded_style_defaults is None:
        nonbonded_style_defaults = 'pair_style lj/cut/coul/long 9.0 9.0\nkspace_style pppm 1e-6\n'

    system, prefix, gro_in, top_in = _load_gromacs(basefile)
    if title is not None:
        system.name = title
    lammps_template_input = os.path.splitext(basefile[0])[0] + ".inp"
    lammps_template_data = os.path.splitext(basefile[0])[0] + ".lmp"
    ll = LammpsParser(lammps_template_input, system=system)
    ll.write(nonbonded_style=nonbonded_style_defaults)

    return lammps_template_input, lammps_template_data

# =============================================================================
def clean_lammps(inp_filename, lmp_filename, gro_filename, top_filename,  ffname):
    """
    This function groups terms in the lmp and inp files.

    Returns:

    """


    lmp_filename_new = os.path.splitext(lmp_filename)+"_clean.lmp"
    inp_filename_new = os.path.splitext(lmp_filename)+"_clean.inp"

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

    energy_terms = _write_lmp_data(lmp_filename_new, ffname, top_filename, lines_top, lines_gro, d_indices_lines)
    _write_lmp_inp(inp_filename_new, ffname, energy_terms)

# =============================================================================
def _write_lmp_data(lmp_filename_new, ffname, top_filename, lines_top, lines_gro, d_indices_lines):

    # Conversion factors
    kjmol_to_kcalmol = 0.239006
    nm_to_angstrom = 10

    natomtypes, atomtypes = _parse_top_atomtypes_section(lines_top, d_indices_lines["atomtypes"], ffname)

    natoms_per_mol, atoms_info =\
        _parse_top_atoms_section(lines_top, d_indices_lines["atoms"])

    nmols = _parse_top_molecules_section(lines_top, d_indices_lines["molecules"])

    nbonds_per_mol, bonds_type, bonds_info =\
        _parse_top_bonds_section(lines_top, d_indices_lines["bonds"], atoms_info)

    nangles_per_mol, angles_type, angles_info =\
        _parse_top_angles_section(lines_top, d_indices_lines["angles"], atoms_info)

    ndihedrals_per_mol, nimpropers_per_mol, dihedrals_type, impropers_type, dihedrals_info = \
        _parse_top_dihedrals_section(lines_top, d_indices_lines["dihedrals"], atoms_info)

    # Check if there is new atom pairs under the section [pairs] of the top file
    _parse_top_pairs_section(lines_top, d_indices_lines["pairs"], dihedrals_info)


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

        # Pair coefficients
        flmp.writelines("\n")
        flmp.writelines("Pair Coeffs\n")
        flmp.writelines("\n")
        for itype in atomtypes:
            # flmp.writelines("{0:<4d} {1:<4d} {2:s} {3:>8.4f} {4:>8.4f}   #{5:s}\n"
            #                 .format(itype.id, itype.id, itype.nbpairs,
            #                         itype.epsilon*kjmol_to_kcalmol,
            #                         itype.sigma*nm_to_angstrom,
            #                         itype.type))
            flmp.writelines("{0:<4d} {1:s} {2:>8.4f} {3:>8.4f}   #{4:s}\n"
                            .format(itype.id, itype.nbpairs,
                                    itype.epsilon*kjmol_to_kcalmol,
                                    itype.sigma*nm_to_angstrom,
                                    itype.type))

        # Bond coefficients
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

        n_energy_terms = {"natoms": natoms_per_mol*nmols,
                          "nbonds": nbonds_per_mol*nmols,
                          "nangles": nangles_per_mol * nmols,
                          "ndihedrals": ndihedrals_per_mol * nmols,
                          "nimpropers": nimpropers_per_mol * nmols}

        return n_energy_terms


# =============================================================================
def _write_lmp_inp(inp_filename_new, ffname, energy_terms, units="real",
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

        if ff.upper() == "OPLSAA" or upper(ff) == "OPLS":
            ffterms = ["lj/cut/coul/long 10.0 10.0", "harmonic", "harmonic", "multi/harmonic", "harmonic"]
            ffterms_special_bonds = ["special_bonds lj 0.0 0.0 0.5 coul 0.0 0.0 0.5"]
            ffterms_pair_modify = ["pair_modify mix geometric tail yes"]
            ffterms_kspace = ["kspace_style pppm 1e-6"]
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
            mass = float(tokens[7])
            ainfo = atoms_info_tuple(idx_atom=idx_atom, type=type, resnumber=resnumber, resname=resname,
                                     atomname=atomname, cgname=cgname, charge=charge, mass=mass)
            atominfo.append(ainfo)
            natoms_per_mol += 1
        idx += 1

    return natoms_per_mol, atominfo

# =============================================================================
def _parse_top_bonds_section(lines_top, idx, atoms_info):

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
                for ak in atoms_info:
                    if ak.idx_atom == idx1:
                        typeA = ak.type
                    elif ak.idx_atom == idx2:
                        typeB = ak.type
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
def _parse_top_angles_section(lines_top, idx, atoms_info):

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
                for ak in atoms_info:
                    if ak.idx_atom == idx1:
                        typeA = ak.type
                    elif ak.idx_atom == idx2:
                        typeB = ak.type
                    elif ak.idx_atom == idx3:
                        typeC = ak.type
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
                print("ERROR: prepare_lammps.py::_parse_top_bonds_section()")
                print("ERROR. [angles] section")
                print("ERROR. Function {} is not yet implemented".format(int(tokens[2])))
                exit()

            nangles += 1
        idx += 1

    return nangles, angletypes, angleinfo

# =============================================================================
def _parse_top_dihedrals_section(lines_top, idx, atoms_info):

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
    while True:
        iline = lines_top[idx]
        if len(iline.replace(" ", "")) < 2:
            break
        if not re.match("^;", iline.strip()) :
            tokens = iline.split()
            # Ryckaert-Bellemans dihedral
            if int(tokens[4]) == 3:
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
                for ak in atoms_info:
                    if ak.idx_atom == idx1:
                        typeA = ak.type
                    elif ak.idx_atom == idx2:
                        typeB = ak.type
                    elif ak.idx_atom == idx3:
                        typeC = ak.type
                    elif ak.idx_atom == idx4:
                        typeD = ak.type
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

            else:
                print("ERROR: prepare_lammps.py::_parse_top_bonds_section()")
                print("ERROR. [dihedrals] section")
                print("ERROR. Function {} is not yet implemented".format(int(tokens[2])))
                exit()
            ndihedrals += 1
        idx += 1

    # for i in dihedraltypes:
    #     print(i)

    return ndihedrals, nimpropers, dihedraltypes, impropertypes, dihedralinfo

# =============================================================================
def _parse_gro_info(lines_gro):

    title = lines_gro[0]
    natoms = int(lines_gro[1])
    coords = np.zeros([natoms, 3])

    # Get coordinates
    for iatom in range(natoms):
        iline = lines_gro[iatom+2]
        coords[iatom, 0] = float(iline[20:29])
        coords[iatom, 1] = float(iline[29:38])
        coords[iatom, 2] = float(iline[38:47])

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

    if ff.upper() == "OPLSAA" or upper(ff) == "OPLS":
        nbpairs = "lj/cut/coul/long"
    else:
        print("ERROR: prepare_lammps.py::_parse_top_atomtypes_section()")
        print("ERROR. [bonds] section")
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


    pass