from collections import defaultdict
import re
import numpy as np

atomic_mass = {'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012,
               'B': 10.811, 'C': 12.011, 'N': 14.007, 'O': 15.999,
               'F': 18.998, 'Ne': 20.180, 'Na': 22.990, 'Mg': 24.305,
               'Al': 26.982, 'Si': 28.086, 'P': 30.974, 'S': 32.066,
               'Cl': 35.453, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078,
               'Sc': 44.956, 'Ti': 47.867, 'V': 50.942, 'Cr': 51.996,
               'Mn': 54.938, 'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693,
               'Cu': 63.546, 'Zn': 65.38, 'Ga': 69.723, 'Ge': 72.631,
               'As': 74.922, 'Se': 78.971, 'Br': 79.904, 'Kr': 84.798,
               'Rb': 84.468, 'Sr': 87.62, 'Y': 88.906, 'Zr': 91.224,
               'Nb': 92.906, 'Mo': 95.95, 'Tc': 98.907, 'Ru': 101.07,
               'Rh': 102.906, 'Pd': 106.42, 'Ag': 107.868, 'Cd': 112.414,
               'In': 114.818, 'Sn': 118.711, 'Sb': 121.760, 'Te': 126.7,
               'I': 126.904, 'Xe': 131.294, 'Cs': 132.905, 'Ba': 137.328,
               'La': 138.905, 'Ce': 140.116, 'Pr': 140.908, 'Nd': 144.243,
               'Pm': 144.913, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25,
               'Tb': 158.925, 'Dy': 162.500, 'Ho': 164.930, 'Er': 167.259,
               'Tm': 168.934, 'Yb': 173.055, 'Lu': 174.967, 'Hf': 178.49,
               'Ta': 180.948, 'W': 183.84, 'Re': 186.207, 'Os': 190.23,
               'Ir': 192.217, 'Pt': 195.085, 'Au': 196.967, 'Hg': 200.592,
               'Tl': 204.383, 'Pb': 207.2, 'Bi': 208.980, 'Po': 208.982,
               'At': 209.987, 'Rn': 222.081, 'Fr': 223.020, 'Ra': 226.025,
               'Ac': 227.028, 'Th': 232.038, 'Pa': 231.036, 'U': 238.029,
               'Np': 237, 'Pu': 244, 'Am': 243, 'Cm': 247, 'Bk': 247,
               'Ct': 251, 'Es': 252, 'Fm': 257, 'Md': 258, 'No': 259,
               'Lr': 262, 'Rf': 261, 'Db': 262, 'Sg': 266, 'Bh': 264,
               'Hs': 269, 'Mt': 268, 'Ds': 271, 'Rg': 272, 'Cn': 285,
               'Nh': 284, 'Fl': 289, 'Mc': 288, 'Lv': 292, 'Ts': 294,
               'Og': 294}
atomic_number = dict(H=1, He=2, Li=3, Be=4, B=5, C=6,
                     N=7, O=8, F=9, Ne=10, Na=11, Mg=12,
                     Al=13, Si=14, P=15, S=16, Cl=17, Ar=18,
                     K=19, Ca=20, Sc=21, Ti=22, V=23, Cr=24,
                     Mn=25, Fe=26, Co=27, Ni=28, Cu=29, Zn=30,
                     Ga=31, Ge=32, As=33, Se=34, Br=35, Kr=36,
                     Rb=37, Sr=38, Y=39, Zr=40, Nb=41, Mo=42,
                     Tc=43, Ru=44, Rh=45, Pd=46, Ag=47,
                     Cd=48, In=49, Sn=50, Sb=51, Te=52,
                     I=53, Xe=54, Cs=55, Ba=56, La=57,
                     Ce=58, Pr=59, Nd=60, Pm=61, Sm=62,
                     Eu=63, Gd=64, Tb=65, Dy=66, Ho=67,
                     Er=68, Tm=69, Yb=70, Lu=71, Hf=72,
                     Ta=73, W=74, Re=75, Os=76, Ir=77,
                     Pt=78, Au=79, Hg=80, Tl=81, Pb=82,
                     Bi=83, Po=84, At=85, Rn=86, Fr=87,
                     Ra=88, Ac=89, Th=90, Pa=91, U=92,
                     Np=93, Pu=94, Am=95, Cm=96, Bk=97,
                     Cf=98, Es=99, Fm=100, Md=101, No=102,
                     Lr=103, Rf=104, Db=105, Sg=106, Bh=107,
                     Hs=108, Mt=109)

atnumber_to_element = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C',
                       7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg',
                       13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar',
                       19: 'K', 20: 'Ca', 21: 'Sc'}


# =============================================================================
def extract_itp_parameters(itp_files_list,
                           start_exclude_list, end_exclude_list, wk="./", log=None):
    """
    Extract all paremeters from itp files
    :param itp_files_list:
    :param start_exclude_list
    :param end_exclude_list
    :param wk
    :param log
    :return:
    """

    # Sort itp_files_list
    # forcefield.itp --> first
    # ffbonded.itp   --> second
    # ffnonbonded    --> third
    # The rest of itps at the end
    tmp = list()
    is_ffbonded_itp = False
    if "./forcefield/forcefield.itp" in itp_files_list:
        itp_files_list.remove("./forcefield/forcefield.itp")
        tmp.append("./forcefield/forcefield.itp")
    if "./forcefield/ffbonded.itp" in itp_files_list:
        itp_files_list.remove("./forcefield/ffbonded.itp")
        tmp.append("./forcefield/ffbonded.itp")
        is_ffbonded_itp = True
    if "./forcefield/ffnonbonded.itp" in itp_files_list:
        itp_files_list.remove("./forcefield/ffnonbonded.itp")
        tmp.append("./forcefield/ffnonbonded.itp")
    itp_files_list = tmp + itp_files_list

    if start_exclude_list is None:
        start_exclude_list = []
    if end_exclude_list is None:
        end_exclude_list = []

    bondtype_dict = defaultdict(list)
    angletype_dict = defaultdict(list)
    dihedraltype_dict = defaultdict(list)
    impropertype_dict = defaultdict(list)
    improper_idx_list = list()
    atomtypes = defaultdict()

    idx_to_atomtype_dict = defaultdict()

    for itpfile in itp_files_list:

        is_atoms_label = False
        is_atomtypes_label = False
        is_bondtypes_label = False
        is_angletypes_label = False
        is_dihedraltypes_label = False
        is_bonds_label = False
        is_angles_label = False
        is_dihedrals_label = False
        is_defaults_label = False

        d_indices_lines = defaultdict()
        with open(itpfile, 'r') as fitp:
            lines_itp = fitp.readlines()
            fitp.seek(0)
            for num, line in enumerate(fitp, 1):
                if re.match(r"^\[(\s*)defaults(\s*)]", line):
                    d_indices_lines["defaults"] = num - 1
                    is_defaults_label = True
                if re.match(r"^\[(\s*)atomtypes(\s*)]", line):
                    d_indices_lines["atomtypes"] = num - 1
                    is_atomtypes_label = True
                if re.match(r"^\[(\s*)bondtypes(\s*)]", line):
                    d_indices_lines["bondtypes"] = num - 1
                    is_bondtypes_label = True
                if re.match(r"^\[(\s*)angletypes(\s*)]", line):
                    d_indices_lines["angletypes"] = num - 1
                    is_angletypes_label = True
                if re.match(r"^\[(\s*)dihedraltypes(\s*)]", line):
                    d_indices_lines["dihedraltypes"] = num - 1
                    is_dihedraltypes_label = True
                if re.match(r"^\[(\s*)nonbond_params(\s*)]", line):
                    d_indices_lines["nonbond_params"] = num - 1
                if re.match(r"^\[(\s*)atoms(\s*)]", line):
                    d_indices_lines["atoms"] = num - 1
                    is_atoms_label = True
                if re.match(r"^\[(\s*)bonds(\s*)]", line):
                    d_indices_lines["bonds"] = num - 1
                    is_bonds_label = True
                if re.match(r"^\[(\s*)angles(\s*)]", line):
                    d_indices_lines["angles"] = num - 1
                    is_angles_label = True
                if re.match(r"^\[(\s*)dihedrals(\s*)]", line):
                    d_indices_lines["dihedrals"] = num - 1
                    is_dihedrals_label = True
                if re.match(r"^\[(\s*)pairs(\s*)]", line):
                    d_indices_lines["pairs"] = num - 1
                if re.match(r"^\[(\s*)moleculetype(\s*)]", line):
                    d_indices_lines["moleculetype"] = num - 1
                if re.match(r"^\[(\s*)molecules(\s*)]", line):
                    d_indices_lines["molecules"] = num - 1
        #
        # # Idx to defaults =================================
        # if is_defaults_label:
        #     idx = d_indices_lines["defaults"] + 1
        #     idx_to_defaults_dict = defaultdict()
        #     while True:
        #         if idx >= len(lines_itp):
        #             break
        #         iline = lines_itp[idx]
        #         if len(iline.replace(" ", "")) < 2:
        #             break
        #         if not re.match("^;", iline.strip()):
        #             tokens = iline.split()
        #             nbfunc = int(tokens[0])
        #             comb_rule = int(tokens[1])
        #             gen_pairs = tokens[2]
        #             fudge_lj = float(tokens[3])
        #             fudge_qq = float(tokens[4])
        #             idx_to_defaults_dict["defaults"] = [nbfunc, comb_rule, gen_pairs, fudge_lj, fudge_qq]
        #         idx += 1

        # # Idx to atomtype =================================
        # if is_atoms_label:
        #     idx = d_indices_lines["atoms"] + 1
        #     idx_to_atomtype_dict = defaultdict()
        #     while True:
        #         if idx >= len(lines_itp):
        #             break
        #         iline = lines_itp[idx]
        #         if len(iline.replace(" ", "")) < 2:
        #             break
        #         if re.match(r'^\[', iline.strip()):
        #             break
        #         if not re.match("^;", iline.strip()):
        #             tokens = iline.split()
        #             iat = int(tokens[0])
        #             itype = tokens[1]
        #             idx_to_atomtype_dict[iat] = itype
        #         idx += 1

        # # Atomtypes  =================================
        # if is_atomtypes_label:
        #     idx = d_indices_lines["atomtypes"] + 1
        #     atomtypes_notfound = list()
        #     while True:
        #         if idx >= len(lines_itp):
        #             break
        #         iline = lines_itp[idx]
        #         if len(iline.replace(" ", "")) < 2:
        #             break
        #         if not re.match("^;", iline.strip()):
        #             tokens = iline.split()
        #             try:
        #                 element = int(tokens[2])
        #                 name = tokens[1]
        #                 charge = float(tokens[4])
        #             except KeyError:
        #                 atomtypes_notfound.append(tokens[0])
        #                 idx += 1
        #                 continue
        #             lines = "{0:<12s}{1:>6s}{2:>3d}{3:>14.5f}{4:>12.5f}{5:>7s}{6:14.5e}{7:14.5e}". \
        #                 format(tokens[0], name, element, float(tokens[3]), charge, tokens[5],
        #                        float(tokens[6]), float(tokens[7]))
        #             atomtypes[tokens[0]] = lines.split()
        #             atomtypes[tokens[0]][2] = int(atomtypes[tokens[0]][2])
        #             atomtypes[tokens[0]][3:5] = [float(item) for item in atomtypes[tokens[0]][3:5]]
        #             atomtypes[tokens[0]][6:] = [float(item) for item in atomtypes[tokens[0]][6:]]
        #
        #         idx += 1

        # # BondTypes =================================
        # if is_bondtypes_label:
        #     idx = d_indices_lines["bondtypes"] + 1
        #     while True:
        #         if idx >= len(lines_itp):
        #             break
        #         iline = lines_itp[idx]
        #         if len(iline.replace(" ", "")) < 2:
        #             break
        #         if not re.match("^;", iline.strip()):
        #             tokens = iline.split()
        #             label = tokens[0]+"-"+tokens[1]
        #             function_bond = int(tokens[2])
        #             r0 = float(tokens[3])
        #             kbond = float(tokens[4])
        #             bondtype_dict[label].append([r0, kbond, function_bond])
        #         idx += 1

        # # AngleTypes =================================
        # if is_angletypes_label:
        #     idx = d_indices_lines["angletypes"] + 1
        #     while True:
        #         if idx >= len(lines_itp):
        #             break
        #         iline = lines_itp[idx]
        #         if len(iline.replace(" ", "")) < 2:
        #             break
        #         if not re.match("^;", iline.strip()):
        #             tokens = iline.split()
        #             itype = tokens[0]
        #             jtype = tokens[1]
        #             ktype = tokens[2]
        #             label = itype + "-" + jtype + "-" + ktype
        #             function_angle = int(tokens[3])
        #             if function_angle == 1 or function_angle == 2:    # Harmonic and cosine, respectively
        #                 theta0 = float(tokens[4])
        #                 kangle = float(tokens[5])
        #                 angletype_dict[label].append([theta0, kangle, function_angle])
        #             elif function_angle == 5:                         # Urey-Bradley
        #                 theta0 = float(tokens[4])
        #                 kangle = float(tokens[5])
        #                 r13 = float(tokens[6])
        #                 kub = float(tokens[7])
        #                 angletype_dict[label].append([theta0, kangle, r13, kub, function_angle])
        #         idx += 1

        # # Dihedral Types =================================
        # if is_dihedraltypes_label:
        #     idx = d_indices_lines["dihedraltypes"] + 1
        #     while True:
        #         if idx >= len(lines_itp):
        #             break
        #         iline = lines_itp[idx]
        #         if len(iline.replace(" ", "")) < 2:
        #             break
        #         if not re.match("^;", iline.strip()):
        #             tokens = iline.split()
        #             itype = int(tokens[0])
        #             jtype = int(tokens[1])
        #             ktype = int(tokens[2])
        #             ltype = int(tokens[3])
        #             labelf = itype + "-" + jtype + "-" + ktype + "-" + ltype
        #             labelr = ltype + "-" + ktype + "-" + jtype + "-" + itype
        #             function_dihedral = int(tokens[4])
        #
        #             if function_dihedral == 9:     # Proper dihedral
        #                 kdih_list = list()
        #                 phase_list = list()
        #                 for j in range(0, 9):
        #                     if len(tokens) < 1:
        #                         break
        #                     kdih_list.append(float(tokens[6]))
        #                     phase_list.append(float(tokens[5]))
        #                     idx += 1
        #                     iline = lines_itp[idx]
        #                     tokens = iline.split()
        #
        #         idx += 1

        # # Bonds ===========================================
        # """
        # Equilibrium distances and force constants can be different due to the QForce fitting procedure
        # """
        # if is_bonds_label:
        #     idx = d_indices_lines["bonds"] + 1
        #     while True:
        #         if idx >= len(lines_itp):
        #             break
        #         iline = lines_itp[idx]
        #         if len(iline.replace(" ", "")) < 2:
        #             break
        #         if not re.match("^;", iline.strip()):
        #             tokens = iline.split()
        #             iat = int(tokens[0])
        #             jat = int(tokens[1])
        #             if iat in start_exclude_list or \
        #                     iat in end_exclude_list or \
        #                     jat in start_exclude_list or \
        #                     jat in end_exclude_list:
        #                 idx += 1
        #                 continue
        #             function_bond = int(tokens[2])
        #             itype = idx_to_atomtype_dict[iat]
        #             jtype = idx_to_atomtype_dict[jat]
        #             labelf = itype + "-" + jtype
        #             labelr = jtype + "-" + itype
        #             try:
        #                 r0 = float(tokens[3])
        #                 kbond = float(tokens[4])
        #                 if labelf not in bondtype_dict and \
        #                         labelr not in bondtype_dict:
        #                     bondtype_dict[labelf].append([r0, kbond, function_bond])
        #                 else:
        #                     if labelf in bondtype_dict:
        #                         lold = labelf
        #                     if labelr in bondtype_dict:
        #                         lold = labelr
        #                     bondtype_dict[lold].append([r0, kbond, function_bond])
        #                 if is_ffbonded_itp:
        #                     m = "\n\t\t WARNING: There are paremeters for {} or {} bond\n".format(labelf, labelr)
        #                     m += "\t\t          in two files: {} and ffbonded.itp\n".format(itpfile)
        #                     m += "\t\t          Please check carefully the parameters.\n".format(itpfile)
        #                     print(m) if log is None else log.warn(m)
        #             except IndexError:
        #                 pass
        #
        #         idx += 1
        #
        # # Angles ===========================================
        # """
        # Equilibrium angles and force constants can be different due to the QForce fitting procedure
        # """
        # if is_angles_label:
        #     idx = d_indices_lines["angles"] + 1
        #     while True:
        #         if idx >= len(lines_itp):
        #             break
        #         iline = lines_itp[idx]
        #         if len(iline.replace(" ", "")) < 2:
        #             break
        #         if not re.match("^;", iline.strip()):
        #             tokens = iline.split()
        #             iat = int(tokens[0])
        #             jat = int(tokens[1])
        #             kat = int(tokens[2])
        #             if iat in start_exclude_list or \
        #                     iat in end_exclude_list or \
        #                     jat in start_exclude_list or \
        #                     jat in end_exclude_list or \
        #                     kat in start_exclude_list or \
        #                     kat in end_exclude_list:
        #                 idx += 1
        #                 continue
        #             function_angle = int(tokens[3])
        #             itype = idx_to_atomtype_dict[iat]
        #             jtype = idx_to_atomtype_dict[jat]
        #             ktype = idx_to_atomtype_dict[kat]
        #             labelf = itype + "-" + jtype + "-" + ktype
        #             labelr = ktype + "-" + jtype + "-" + itype
        #             try:
        #                 if function_angle == 1 or function_angle == 2:    # Harmonic and cosine forms
        #                     theta0 = float(tokens[4])
        #                     kangle = float(tokens[5])
        #                     if labelf not in angletype_dict and \
        #                             labelr not in angletype_dict:
        #                         angletype_dict[labelf].append([theta0, kangle, function_angle])
        #                     else:
        #                         if labelf in angletype_dict:
        #                             lold = labelf
        #                         if labelr in angletype_dict:
        #                             lold = labelr
        #                         angletype_dict[lold].append([theta0, kangle, function_angle])
        #                 elif function_angle == 5:                         # Urey-Bradley
        #                     theta0 = float(tokens[4])
        #                     kangle = float(tokens[5])
        #                     r13 = float(tokens[6])
        #                     kub = float(tokens[7])
        #                     if labelf not in angletype_dict and \
        #                             labelr not in angletype_dict:
        #                         angletype_dict[labelf].append([theta0, kangle, r13, kub, function_angle])
        #                     else:
        #                         if labelf in angletype_dict:
        #                             lold = labelf
        #                         if labelr in angletype_dict:
        #                             lold = labelr
        #                         angletype_dict[lold].append([theta0, kangle, r13, kub, function_angle])
        #                 else:
        #                     m = "\n\t\t\t ERROR Angle function {} does not exist".format(function_angle)
        #                     print(m) if log is None else log.error(m)
        #                     exit()
        #                 if is_ffbonded_itp:
        #                     m = "\n\t\t WARNING: There are paremeters for {} or {} angle\n".format(labelf, labelr)
        #                     m += "\t\t          in two files: {} and ffbonded.itp\n".format(itpfile)
        #                     m += "\t\t          Please check carefully the parameters.\n".format(itpfile)
        #                     print(m) if log is None else log.warn(m)
        #             except IndexError:
        #                 pass
        #
        #         idx += 1

        # Dihedrals ===========================================
        """
        Equilibrium dihedrals and force constants can be different due to the QForce fitting procedure
        """
        if is_dihedrals_label:
            idx = d_indices_lines["dihedrals"] + 1
            while True:
                if idx >= len(lines_itp):
                    break
                iline = lines_itp[idx]
                if len(iline.replace(" ", "")) < 2:
                    break
                if not re.match("^;", iline.strip()):
                    tokens = iline.split()
                    iat = int(tokens[0])
                    jat = int(tokens[1])
                    kat = int(tokens[2])
                    lat = int(tokens[3])
                    if iat in start_exclude_list or \
                            iat in end_exclude_list or \
                            jat in start_exclude_list or \
                            jat in end_exclude_list or \
                            kat in start_exclude_list or \
                            kat in end_exclude_list or \
                            lat in start_exclude_list or \
                            lat in end_exclude_list:
                        idx += 1
                        continue
                    itype = idx_to_atomtype_dict[iat]
                    jtype = idx_to_atomtype_dict[jat]
                    ktype = idx_to_atomtype_dict[kat]
                    ltype = idx_to_atomtype_dict[lat]
                    labelf = itype + "-" + jtype + "-" + ktype + "-" + ltype
                    labelr = ltype + "-" + ktype + "-" + jtype + "-" + itype

                    function_dihedral = int(tokens[4])
                    try:
                        if function_dihedral == 9:     # Proper dihedral
                            kdih_list = list()
                            phase_list = list()
                            for j in range(0, 9):
                                if len(tokens) < 1:
                                    break
                                kdih_list.append(float(tokens[6]))
                                phase_list.append(float(tokens[5]))
                                idx += 1
                                iline = lines_itp[idx]
                                tokens = iline.split()
                        elif function_dihedral == 3:    # Proper Ryckaert-Bellemans dihedral
                            kdih_list = [float(i) for i in tokens[5:]]
                            idx += 1
                        elif function_dihedral == 2:    # Improper Harmonic type
                            impr0 = float(tokens[5])
                            kimp = float(tokens[6])
                            impr_param = [impr0, kimp, 0.0, function_dihedral]
                            idx += 1
                        elif function_dihedral == 4:   # Improper periodic type
                            impr0 = float(tokens[5])
                            kimp = float(tokens[6])
                            multi = int(tokens[7])
                            impr_param = [impr0, kimp, multi, function_dihedral]
                            idx += 1
                        else:
                            m = "\n\t\t\t ERROR Dihedral function {} does not exist".format(function_dihedral)
                            print(m) if log is None else log.error(m)
                            exit()
                        if function_dihedral != 2 and function_dihedral != 4:     # Only proper dihedrals
                            if labelf not in dihedraltype_dict and \
                                    labelr not in dihedraltype_dict:
                                dihedraltype_dict[labelf].append([kdih_list, function_dihedral])
                            else:
                                if labelf in dihedraltype_dict:
                                    lold = labelf
                                if labelr in dihedraltype_dict:
                                    lold = labelr
                                dihedraltype_dict[lold].append([kdih_list, function_dihedral])
                        else:
                            if labelf not in impropertype_dict and \
                                    labelr not in impropertype_dict:
                                impropertype_dict[labelf].append(impr_param)
                            else:
                                if labelf in impropertype_dict:
                                    lold = labelf
                                if labelr in impropertype_dict:
                                    lold = labelr
                                impropertype_dict[lold].append(impr_param)
                            improper_idx_list.append([iat, jat, kat, lat])
                    except IndexError:
                        idx += 1
                        pass
                else:
                    idx += 1

    m = "\t\t  ===== ATOM TYPES =====\n"
    for key, item in atomtypes.items():
        m += "\t\t Atom type {0:14s} ({1:3s}): sigma = {2:7.4f} +- 0.000 nm (0.00 %) " \
             ";epsilon = {3:8.4f} +- 0.0 kJ/mol (0.0%) \n". \
            format(key, item[1], item[6], item[7])
    atomtypes_notfound = set(atomtypes_notfound)
    for item in atomtypes_notfound:
        m += "\t\t WARNING: Atom type {} info is not available in the info ff file.\n".format(item)
    print(m) if log is None else log.info(m)

    # Averages ==================================================
#cJ    bondtype_avg_dict = defaultdict()
    for key, coeff_list in bondtype_dict.items():
        r0_list = []
        kbond_list = []
        function_list = []
        for item in coeff_list:
            r0_list.append(item[0])
            kbond_list.append(item[1])
            function_list.append(item[2])
        m1 = np.mean(r0_list)
        std1 = np.std(r0_list)
        m2 = np.mean(kbond_list)
        std2 = np.std(kbond_list)

        bondtype_avg_dict[key] = [m1, std1, m2, std2, len(r0_list), function_list[-1]]

#cJ    angletype_avg_dict = defaultdict()
    for key, coeff_list in angletype_dict.items():
        theta0_list = []
        kangle_list = []
        r13_list = []
        kub_list = []
        if coeff_list[0][2] == 1 or coeff_list[0][2] == 2:
            for item in coeff_list:
                theta0_list.append(item[0])
                kangle_list.append(item[1])
                function_list.append(item[2])
            m1 = np.mean(theta0_list)
            std1 = np.std(theta0_list)
            m2 = np.mean(kangle_list)
            std2 = np.std(kangle_list)
            angletype_avg_dict[key] = [m1, std1, m2, std2, len(theta0_list), function_list[-1]]
        elif coeff_list[0][4] == 5:
            for item in coeff_list:
                theta0_list.append(item[0])
                kangle_list.append(item[1])
                r13_list.append(item[2])
                kub_list.append(item[3])
                function_list.append(item[4])
            m1 = np.mean(theta0_list)
            std1 = np.std(theta0_list)
            m2 = np.mean(kangle_list)
            std2 = np.std(kangle_list)
            m3 = np.mean(r13_list)
            std3 = np.std(r13_list)
            m4 = np.mean(kub_list)
            std4 = np.std(kub_list)
            angletype_avg_dict[key] = [m1, std1, m2, std2, m3, std3, m4, std4, len(theta0_list), function_list[-1]]

#cJ    dihedraltype_avg_dict = defaultdict()
#cJ    dihedraltype_std_dict = defaultdict()
#cJ    dihedraltype_function_dict = defaultdict()
    for key, coeff_list in dihedraltype_dict.items():
        ncoeff = len(coeff_list[0][0])
        nsamples = 0
        acc_coeffs = [0.0 for _ in range(0, ncoeff)]
        std_coeffs = [0.0 for _ in range(0, ncoeff)]
        for icoeff in coeff_list:
            ncoeff_new = len(icoeff[0])
            # Check that all samples for the same dihedral have the same number of coefficients
            if ncoeff_new != ncoeff:
                m = "\n\t\tERROR: The number of coefficients must be the same for each sample!!!\n"
                m += "\t\tERROR: Dihedral problem {}\n".format(key)
                m += "\t\tERROR:{}".format([i for i in coeff_list])
                m += "\n"
                print(m) if log is None else log.info(m)
                exit()
            # Accumulate coefficients
            for jj in range(0, len(icoeff[0])):
                acc_coeffs[jj] += icoeff[0][jj]
            nsamples += 1

        dihedraltype_avg_dict[key] = [acc_coeffs[j] / nsamples for j in range(0, len(acc_coeffs))]
        dihedraltype_function_dict[key] = icoeff[1]

        for icoeff in coeff_list:
            # Accumulate coefficients for std
            for jj in range(0, len(icoeff)):
                std_coeffs[jj] += (icoeff[0][jj] - dihedraltype_avg_dict[key][jj]) ** 2
            nsamples += 1

        dihedraltype_std_dict[key] = [np.sqrt(std_coeffs[j] / (nsamples - 1)) for j in range(0, len(std_coeffs))]

#cJ   impropertype_avg_dict = defaultdict()
    for key, coeff_list in impropertype_dict.items():
        impr0_list = []
        kimp_list = []
        multiplicity_list = []

        for item in coeff_list:
            impr0_list.append(item[0])
            kimp_list.append(item[1])
            multiplicity_list.append(item[2])
            function_list.append(item[-1])
        m1 = np.mean(impr0_list)
        std1 = np.std(impr0_list)
        m2 = np.mean(kimp_list)
        std2 = np.std(kimp_list)
        m3 = np.mean(multiplicity_list)
        std3 = np.std(multiplicity_list)
        impropertype_avg_dict[key] = [m1, std1, m2, std2, m3, std3, len(impr0_list), function_list[-1]]

    m = "\t\t  ===== BOND TYPES =====\n"
    for key, coeff_list in bondtype_avg_dict.items():
        p1 = coeff_list[1] / coeff_list[0] * 100.
        p2 = coeff_list[3] / coeff_list[2] * 100.
        m += "\t\t Bond type {0:30s}: r0 = {1:7.4f} +- {2:7.4f} nm ({3:5.2f} %) " \
             ";k_b = {4:9.1f} +- {5:9.1f} kJ/nm^2 ({6:5.2f}%) ({7:d} elements)\n". \
            format(key, coeff_list[0], coeff_list[1], p1, coeff_list[2], coeff_list[3], p2, coeff_list[4])

    print(m) if log is None else log.info(m)

    m = "\t\t  ===== ANGLE TYPES =====\n"
    for key, coeff_list in angletype_avg_dict.items():
        p1 = coeff_list[1] / coeff_list[0] * 100.
        p2 = coeff_list[3] / coeff_list[2] * 100.
        function_angle = coeff_list[-1]
        if function_angle == 1 or function_angle == 2:   # Harmonic or Cosine forms
            m += "\t\t Angle type {0:35s}: theta0 = {1:<7.1f} +- {2:7.1f} deg ({3:5.2f} %) " \
                 ";k_b = {4:9.1f} +- {5:9.1f} kJ/rad^2 ({6:5.2f}%) ({7:d} elements)\n". \
                format(key, coeff_list[0], coeff_list[1], p1, coeff_list[2], coeff_list[3], p2, coeff_list[4])
        elif function_angle == 5:
            if coeff_list[4] != 0.0:
                p3 = coeff_list[5] / coeff_list[4] * 100.
            else:
                p3 = 0.0
            if coeff_list[6] != 0.0:
                p4 = coeff_list[7] / coeff_list[6] * 100.
            else:
                p4 = 0.0
            m += "\t\t Angle type {0:35s}: theta0 = {1:<7.1f} +- {2:7.1f} deg ({3:5.2f} %) " \
                 ";k_b = {4:9.1f} +- {5:9.1f} kJ/rad^2 ({6:5.2f}%) " \
                 ";r13 = {7:9.2f} +- {8:9.2f} nm ({9:5.2f}%)" \
                 ";kub = {10:9.1f} +- {11:9.1f} nm ({12:5.2f}%)" \
                 "  ({13:d} elements)\n". \
                 format(key, coeff_list[0], coeff_list[1], p1,
                        coeff_list[2], coeff_list[3], p2,
                        coeff_list[4], coeff_list[5], p3,
                        coeff_list[6], coeff_list[7], p4,
                        coeff_list[8])

    print(m) if log is None else log.info(m)

    m = "\t\t  ===== DIHEDRAL TYPES =====\n"
    for key, coeff_list in dihedraltype_avg_dict.items():
        m += "\t\t Dihedral type {0:45s}: ". \
            format(key)
        for jdx, i in enumerate(coeff_list):
            m += "{0:7.3f} +- {1:7.3f} ;".format(i, dihedraltype_std_dict[key][jdx])
        m += "\n"
    print(m) if log is None else log.info(m)

    m = "\t\t  ===== IMPROPER TYPES =====\n"
    for key, coeff_list in impropertype_avg_dict.items():
        p1 = coeff_list[1] / coeff_list[0] * 100.
        p2 = coeff_list[3] / coeff_list[2] * 100.
        function = coeff_list[-1]
        if function == 2:
            m += "\t\t Harmonic Improprer type {0:45s}: impr0 = {1:<7.1f} +- {2:7.1f} deg ({3:5.2f} %) " \
                 ";k_impr = {4:9.1f} +- {5:9.1f} kJ/rad^2 ({6:5.2f}%) ({7:d} elements)\n". \
                format(key, coeff_list[0], coeff_list[1], p1, coeff_list[2], coeff_list[3], p2, coeff_list[-1])
        elif function == 4:
            m += "\t\t Periodic Improprer type {0:45s}: impr0 = {1:<7.1f} +- {2:7.1f} deg ({3:5.2f} %) " \
                 ";k_impr = {4:9.1f} +- {5:9.1f} kJ ({6:5.2f}%) ({7:d} elements)\n". \
                format(key, coeff_list[0], coeff_list[1], p1, coeff_list[2], coeff_list[3], p2, coeff_list[-1])

    print(m) if log is None else log.info(m)

    return idx_to_defaults_dict, bondtype_avg_dict, \
        angletype_avg_dict, dihedraltype_avg_dict, dihedraltype_function_dict, impropertype_avg_dict, \
        improper_idx_list, atomtypes


# =============================================================================
#def extract_gromacs_parameters()