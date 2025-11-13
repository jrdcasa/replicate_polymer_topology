from collections import defaultdict
import re

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
def extract_gromacs_parameters(itp_files_list,
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
    forcefield_itp_list = list()
    is_ffbonded_itp = False
    if "./forcefield/forcefield.itp" in itp_files_list:
        itp_files_list.remove("./forcefield/forcefield.itp")
        forcefield_itp_list.append("./forcefield/forcefield.itp")
    if "./forcefield/ffbonded.itp" in itp_files_list:
        itp_files_list.remove("./forcefield/ffbonded.itp")
        forcefield_itp_list.append("./forcefield/ffbonded.itp")
        is_ffbonded_itp = True
    if "./forcefield/ffnonbonded.itp" in itp_files_list:
        itp_files_list.remove("./forcefield/ffnonbonded.itp")
        forcefield_itp_list.append("./forcefield/ffnonbonded.itp")

    # Setup dictionaries
    improper_idx_list = list()

    bondtype_avg_dict = defaultdict()
    angletype_avg_dict = defaultdict()
    dihedraltype_avg_dict = defaultdict()
    dihedraltype_std_dict = defaultdict()
    dihedraltype_function_dict = defaultdict()
    impropertype_avg_dict = defaultdict()
    idx_to_atomtype_dict = defaultdict()

    d_indices_lines = defaultdict(list)

    # [ defaults ] label should be in the forcefield.itp file ===================================
    ff_name = "./forcefield/forcefield.itp"
    try:
        idx_to_defaults_dict = aux_get_idx_to_defaults_dict(ff_name, d_indices_lines, log)
    except FileNotFoundError as e:
        m = "\n\t\tERROR file {} does not exist.\n".format(ff_name)
        m += "\t\t{}".format(e)
        print(m) if log is None else log.error(m)
        exit()

    # [ atomtypes ] label should be in the ffnonbonded.itp file ==================================
    ff_name = "./forcefield/ffnonbonded.itp"
    try:
        atomtypes = aux_get_atomtypes(ff_name, d_indices_lines, log)
    except FileNotFoundError as e:
        m = "\n\t\tERROR file {} does not exist.\n".format(ff_name)
        m += "\t\t{}".format(e)
        print(m) if log is None else log.error(m)
        exit()

    # [ bondtypes ] label should be in the ffbonded.itp file ==================================
    # [ angletypes ] label should be in the ffbonded.itp file =================================
    # [ dihedraltypes ] label should be in the ffbonded.itp file ==============================
    # [ pairtypes ] label should be in the ffbonded.itp file ==================================
    ff_name = "./forcefield/ffbonded.itp"
    try:
        bondtype_dict, angletype_dict, dihedraltype_dict, \
            impropertype_dict, pairtype_dict = aux_get_bondedtypes(ff_name, d_indices_lines, log)
    except FileNotFoundError as e:
        bondtype_dict = defaultdict(list)
        angletype_dict = defaultdict(list)
        dihedraltype_dict = defaultdict(list)
        impropertype_dict = defaultdict(list)
        pairtype_dict = defaultdict(list)

    for itpfile in itp_files_list:

        improper_idx_list, bondtype_repeat_dict, angletype_repeat_dict, \
            dihedraltype_repeat_dict, impropertype_repeat_dict, pairtype_repeat_dict =\
            aux_get_moleculetype(itpfile, d_indices_lines, bondtype_dict,
                                 angletype_dict, dihedraltype_dict, impropertype_dict, pairtype_dict, log)

    m = "\t\t  ===== BOND TYPES =====\n"
    for key, coeff_list in bondtype_dict.items():
        m += "\t\t Bond type {0:30s}: function={1:2d} coeffs=".format(key, coeff_list[0][-1])
        for icoeff in coeff_list[0][:-1]:
            m += "{0:12.4f}, ".format(icoeff)
        m = m[:-2] + "\n"

    print(m) if log is None else log.info(m)

    m = "\t\t  ===== ANGLE TYPES =====\n"
    for key, coeff_list in angletype_dict.items():
        m += "\t\t Angle type {0:40s}: function={1:2d} coeffs=".format(key, coeff_list[0][-1])
        for icoeff in coeff_list[0][:-1]:
            m += "{0:12.4f}, ".format(icoeff)
        m = m[:-2] + "\n"

    print(m) if log is None else log.info(m)

    m = "\t\t  ===== DIHEDRAL TYPES =====\n"
    for key, coeff_list in dihedraltype_dict.items():
        #  Multiline Dihedrals
        if coeff_list[0][-1] == 9:
            for idx in range(len(coeff_list)):
                m += "\t\t Dihedral type {0:50s}: function={1:2d} coeffs=".format(key, coeff_list[0][-1])
                for jcoeff in coeff_list[idx][0:-2]:
                    m += "{0:12.4f}, ".format(jcoeff)
                m += "{0:2d}".format(coeff_list[idx][-2])
                m += "\n"
        else:
            m += "\t\t Dihedral type {0:50s}: function={1:2d} coeffs=".format(key, coeff_list[0][-1])
            for icoeff in coeff_list[0][:-1]:
                m += "{0:12.4f}, ".format(icoeff)
            m = m[:-2] + "\n"
        dihedraltype_function_dict[key] = coeff_list[0][-1]
    print(m) if log is None else log.info(m)

    m = "\t\t  ===== IMPROPER TYPES =====\n"
    for key, coeff_list in impropertype_dict.items():
        m += "\t\t Improper type {0:50s}: function={1:2d} coeffs=".format(key, coeff_list[0][-1])
        for icoeff in coeff_list[0][:-1]:
            m += "{0:12.4f}, ".format(icoeff)
        m = m[:-2] + "\n"

    m = "\t\t  ===== PAIR TYPES =====\n"
    for key, coeff_list in pairtype_dict.items():
        m += "\t\t Pair type {0:30s}: function={1:2d} coeffs=".format(key, coeff_list[0][-1])
        for icoeff in coeff_list[0][:-1]:
            m += "{0:12.4f}, ".format(icoeff)
        m = m[:-2] + "\n"

    print(m) if log is None else log.info(m)

    return idx_to_defaults_dict, bondtype_dict, \
        angletype_dict, dihedraltype_dict, impropertype_dict, \
        improper_idx_list, atomtypes, dihedraltype_function_dict


# ---------------------------------------------------------------------------------------
def aux_get_idx_to_defaults_dict(ff_name, d_indices_lines, log):
    idx_to_defaults_dict = defaultdict()

    # [ defaults ] label should be in the forcefield.itp file ===================================
    is_defaults_label = False
    with open(ff_name, 'r') as ff_itp:
        ff_itp.seek(0)
        for num, line in enumerate(ff_itp, 1):
            if re.match(r"^\[(\s*)defaults(\s*)]", line):
                d_indices_lines["defaults"].append({"file": ff_name, "idx": num - 1})
                is_defaults_label = True
        if len(d_indices_lines["defaults"]) != 1 or not is_defaults_label:
            m = "\n\t\t ERROR: [ defaults ] label must appear only once in {}\n".format(ff_name)
            m += "\t\t        It appears {} times.\n".format(len(d_indices_lines["defaults"]))
            m += "\t\t        Check the [ defaults ] directives in the input forcefield files.\n"
            print(m) if log is None else log.error(m)
            exit()
        else:
            ff_itp.seek(0)
            lines_itp = ff_itp.readlines()
            idx = d_indices_lines["defaults"][0]["idx"] + 1
            while True:
                if idx >= len(lines_itp):
                    break
                iline = lines_itp[idx]
                if len(iline.replace(" ", "")) < 2:
                    break
                if not re.match("^;", iline.strip()):
                    tokens = iline.split()
                    nbfunc = int(tokens[0])
                    comb_rule = int(tokens[1])
                    gen_pairs = tokens[2]
                    fudge_lj = float(tokens[3])
                    fudge_qq = float(tokens[4])
                    idx_to_defaults_dict["defaults"] = [nbfunc, comb_rule, gen_pairs, fudge_lj, fudge_qq]
                idx += 1

    return idx_to_defaults_dict


# ---------------------------------------------------------------------------------------
def aux_get_atomtypes(ff_name, d_indices_lines, log):
    atomtypes = defaultdict(list)

    # [ atomtypes ] label should be in the ffnonbonded.itp file ==================================
    is_atomtype_label = False
    with open(ff_name, 'r') as ff_itp:
        ff_itp.seek(0)
        for num, line in enumerate(ff_itp, 1):
            if re.match(r"^\[(\s*)atomtypes(\s*)]", line):
                d_indices_lines["atomtypes"].append({"file": ff_name, "idx": num - 1})
                is_atomtype_label = True
        if not is_atomtype_label:
            m = "\n\t\t ERROR: [ atomtypes ] label must appear at least once in {}\n".format(ff_name)
            m += "\t\t        It appears {} times.\n".format(len(d_indices_lines["atomtypes"]))
            m += "\t\t        Check the [ atomtypes ] directives in the input forcefield files.\n"
            print(m) if log is None else log.error(m)
            exit()
        else:
            ff_itp.seek(0)
            lines_itp = ff_itp.readlines()
            idx = d_indices_lines["atomtypes"][0]["idx"] + 1
            atomtypes_notfound = list()
            while True:
                if idx >= len(lines_itp):
                    break
                iline = lines_itp[idx]
                if len(iline.replace(" ", "")) < 2:
                    break
                if not re.match("^;", iline.strip()):
                    tokens = iline.split()
                    try:
                        element = int(tokens[2])
                        bondtypename = tokens[1]
                        charge = float(tokens[4])
                    except KeyError:
                        atomtypes_notfound.append(tokens[0])
                        idx += 1
                        continue
                    lines = "{0:<12s}{1:>6s}{2:>3d}{3:>14.5f}{4:>12.5f}{5:>7s}{6:14.5e}{7:14.5e}". \
                        format(tokens[0], bondtypename, element, float(tokens[3]), charge, tokens[5],
                               float(tokens[6]), float(tokens[7]))
                    if tokens[0] in atomtypes.keys():
                        m = " WARNING: The atomtype {} appears more than once in [ atomtypes ]".format(tokens[0])
                        print(m) if log is None else log.warn(m)
                    atomtypes[tokens[0]] = lines.split()
                    atomtypes[tokens[0]][2] = int(atomtypes[tokens[0]][2])
                    atomtypes[tokens[0]][3:5] = [float(item) for item in atomtypes[tokens[0]][3:5]]
                    atomtypes[tokens[0]][6:] = [float(item) for item in atomtypes[tokens[0]][6:]]

                idx += 1

    return atomtypes


# ---------------------------------------------------------------------------------------
def aux_get_bondedtypes(ff_name, d_indices_lines, log):

    bondtype_dict = defaultdict(list)
    angletype_dict = defaultdict(list)
    dihedraltype_dict = defaultdict(list)
    impropertype_dict = defaultdict(list)
    pairtype_dict = defaultdict(list)

    # [ bondtypes ] label should be in the ffbonded.itp file ==================================
    is_bondtype_label = False
    with open(ff_name, 'r') as ff_itp:
        ff_itp.seek(0)
        for num, line in enumerate(ff_itp, 1):
            if re.match(r"^\[(\s*)bondtypes(\s*)]", line):
                d_indices_lines["bondtypes"].append({"file": ff_name, "idx": num - 1})
                is_bondtype_label = True
        if not is_bondtype_label:
            pass
        else:
            ff_itp.seek(0)
            lines_itp = ff_itp.readlines()
            idx = d_indices_lines["bondtypes"][0]["idx"] + 1
            while True:
                if idx >= len(lines_itp):
                    break
                iline = lines_itp[idx]
                if len(iline.replace(" ", "")) < 2:
                    break
                if not re.match("^;", iline.strip()):
                    tokens = iline.split()
                    label = tokens[0] + "-" + tokens[1]
                    function_bond = int(tokens[2])
                    # r0, kbbond
                    list_parameters = [float(i) for i in tokens[3:]]
                    list_parameters.append(function_bond)
                    bondtype_dict[label].append(list_parameters)
                idx += 1

    # [ angletypes ] label should be in the ffbonded.itp file ==================================
    ff_name = "./forcefield/ffbonded.itp"
    is_angletype_label = False
    with open(ff_name, 'r') as ff_itp:
        ff_itp.seek(0)
        for num, line in enumerate(ff_itp, 1):
            if re.match(r"^\[(\s*)angletypes(\s*)]", line):
                d_indices_lines["angletypes"].append({"file": ff_name, "idx": num - 1})
                is_angletype_label = True
        if not is_angletype_label:
            pass
        else:
            ff_itp.seek(0)
            lines_itp = ff_itp.readlines()
            idx = d_indices_lines["angletypes"][0]["idx"] + 1
            while True:
                if idx >= len(lines_itp):
                    break
                iline = lines_itp[idx]
                if len(iline.replace(" ", "")) < 2:
                    break
                if not re.match("^;", iline.strip()):
                    tokens = iline.split()
                    itype = tokens[0]
                    jtype = tokens[1]
                    ktype = tokens[2]
                    label = itype + "-" + jtype + "-" + ktype
                    function_angle = int(tokens[3])
                    if function_angle == 1 or function_angle == 2:  # Harmonic and cosine, respectively
                        theta0 = float(tokens[4])
                        kangle = float(tokens[5])
                        angletype_dict[label].append([theta0, kangle, function_angle])
                    elif function_angle == 5:  # Urey-Bradley
                        theta0 = float(tokens[4])
                        kangle = float(tokens[5])
                        r13 = float(tokens[6])
                        kub = float(tokens[7])
                        angletype_dict[label].append([theta0, kangle, r13, kub, function_angle])
                idx += 1

    # [ dihedraltypes ] label should be in the ffbonded.itp file ==================================
    ff_name = "./forcefield/ffbonded.itp"
    is_dihedraltype_label = False
    with open(ff_name, 'r') as ff_itp:
        ff_itp.seek(0)
        for num, line in enumerate(ff_itp, 1):
            if re.match(r"^\[(\s*)dihedraltypes(\s*)]", line):
                d_indices_lines["dihedraltypes"].append({"file": ff_name, "idx": num - 1})
                is_dihedraltype_label = True
        if not is_dihedraltype_label:
            pass
        else:
            ff_itp.seek(0)
            lines_itp = ff_itp.readlines()
            idx = d_indices_lines["dihedraltypes"][0]["idx"] + 1
            while True:
                if idx >= len(lines_itp):
                    break
                iline = lines_itp[idx]
                if len(iline.replace(" ", "")) < 2:
                    break
                if not re.match("^;", iline.strip()):
                    tokens = iline.split()
                    itype = tokens[0]
                    jtype = tokens[1]
                    ktype = tokens[2]
                    ltype = tokens[3]
                    labelf = itype + "-" + jtype + "-" + ktype + "-" + ltype
                    labelr = ltype + "-" + ktype + "-" + jtype + "-" + itype
                    function_dihedral = int(tokens[4])
                    if function_dihedral == 9:  # Proper dihedral
                        tokens = iline.split()
                        tokens[5] = float(tokens[5])
                        tokens[6] = float(tokens[6])
                        tokens[7] = int(tokens[7])
                        if labelf in dihedraltype_dict:
                            list_parameters = tokens[5:]
                            list_parameters.append(function_dihedral)
                            dihedraltype_dict[labelf].append(list_parameters)
                        elif labelr in dihedraltype_dict:
                            list_parameters = tokens[5:]
                            list_parameters.append(function_dihedral)
                            dihedraltype_dict[labelr].append(list_parameters)
                        else:
                            list_parameters = tokens[5:]
                            list_parameters.append(function_dihedral)
                            dihedraltype_dict[labelf].append(list_parameters)
                    elif function_dihedral == 3:  # Proper Ryckaert-Bellemans dihedral
                        kdih_list = [float(i) for i in tokens[5:]]
                    elif function_dihedral == 2:  # Improper Harmonic type
                        impr0 = float(tokens[5])
                        kimp = float(tokens[6])
                        impr_param = [impr0, kimp, 0.0, function_dihedral]
                    elif function_dihedral == 4:  # Improper periodic type
                        impr0 = float(tokens[5])
                        kimp = float(tokens[6])
                        multi = int(tokens[7])
                        impr_param = [impr0, kimp, multi, function_dihedral]
                    else:
                        m = "\n\t\t\t ERROR Dihedral function {} does not exist".format(function_dihedral)
                        print(m) if log is None else log.error(m)
                        exit()
                    if function_dihedral != 2 and function_dihedral != 4 and \
                            function_dihedral != 9:  # Only proper dihedrals and not multiline torsion
                        if labelf not in dihedraltype_dict and \
                                labelr not in dihedraltype_dict:
                            dihedraltype_dict[labelf].append([kdih_list, function_dihedral])
                        else:
                            if labelf in dihedraltype_dict:
                                lold = labelf
                            if labelr in dihedraltype_dict:
                                lold = labelr
                            dihedraltype_dict[lold].append([kdih_list, function_dihedral])
                    elif function_dihedral == 9:
                        # The dihedral dict is already setup
                        pass
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
                idx += 1

    # [ pairtypes ] label should be in the ffbonded.itp file ==================================
    ff_name = "./forcefield/ffbonded.itp"
    is_pairtype_label = False
    with open(ff_name, 'r') as ff_itp:
        ff_itp.seek(0)
        for num, line in enumerate(ff_itp, 1):
            if re.match(r"^\[(\s*)pairtypes(\s*)]", line):
                d_indices_lines["pairtypes"].append({"file": ff_name, "idx": num - 1})
                is_pairtype_label = True
        if not is_pairtype_label:
            pass
        else:
            ff_itp.seek(0)
            lines_itp = ff_itp.readlines()
            idx = d_indices_lines["pairtypes"][0]["idx"] + 1
            while True:
                if idx >= len(lines_itp):
                    break
                iline = lines_itp[idx]
                if len(iline.replace(" ", "")) < 2:
                    break
                if not re.match("^;", iline.strip()):
                    tokens = iline.split()
                    label = tokens[0] + "-" + tokens[1]
                    function_pair = int(tokens[2])
                    if function_pair == 2:
                        print("ERROR: Extra LJ or Coulomb potential is not implemented.")
                        exit()
                    sigma14 = float(tokens[3])
                    epsilon14 = float(tokens[4])
                    pairtype_dict[label].append([sigma14, epsilon14, function_pair])
                idx += 1

    return bondtype_dict, angletype_dict, dihedraltype_dict, impropertype_dict, pairtype_dict


# ---------------------------------------------------------------------------------------
def aux_get_moleculetype(ff_name, d_indices_lines, bondtype_dict,
                         angletype_dict, dihedraltype_dict,
                         impropertype_dict, pairtype_dict, log):

    idx_to_atomtype_dict = defaultdict()

    bondtype_repeat_dict = defaultdict(list)
    angletype_repeat_dict = defaultdict(list)
    dihedraltype_repeat_dict = defaultdict(list)
    impropertype_repeat_dict = defaultdict(list)
    pairtype_repeat_dict = defaultdict(list)
    improper_idx_list = list()

    # [ atoms ] label =============================================================
    is_atoms_label = False
    with open(ff_name, 'r') as ff_itp:
        ff_itp.seek(0)
        for num, line in enumerate(ff_itp, 1):
            if re.match(r"^\[(\s*)atoms(\s*)]", line):
                d_indices_lines["atoms"].append({"file": ff_name, "idx": num - 1})
                is_atoms_label = True
        # Idx to atomtype =================================
        if is_atoms_label:
            ff_itp.seek(0)
            lines_itp = ff_itp.readlines()
            idx = d_indices_lines["atoms"][0]["idx"] + 1
            while True:
                if idx >= len(lines_itp):
                    break
                iline = lines_itp[idx]
                if len(iline.replace(" ", "")) < 2:
                    break
                if re.match(r'^\[', iline.strip()):
                    break
                if not re.match("^;", iline.strip()):
                    tokens = iline.split()
                    iat = int(tokens[0])
                    itype = tokens[1]
                    idx_to_atomtype_dict[iat] = itype
                idx += 1

    # [ bonds ] label =============================================================
    is_bonds_label = False
    with open(ff_name, 'r') as ff_itp:
        ff_itp.seek(0)
        for num, line in enumerate(ff_itp, 1):
            if re.match(r"^\[(\s*)bonds(\s*)]", line):
                d_indices_lines["bonds"].append({"file": ff_name, "idx": num - 1})
                is_bonds_label = True
        if is_bonds_label:
            ff_itp.seek(0)
            lines_itp = ff_itp.readlines()
            idx = d_indices_lines["bonds"][0]["idx"] + 1
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
                    function_bond = int(tokens[2])
                    itype = idx_to_atomtype_dict[iat]
                    jtype = idx_to_atomtype_dict[jat]
                    labelf = itype + "-" + jtype
                    labelr = jtype + "-" + itype
                    # The label is not in the bondtype_dict, adding the parameters
                    if labelf not in bondtype_dict and \
                            labelr not in bondtype_dict:
                        # r0, kbbond
                        list_parameters = [float(i) for i in tokens[3:]]
                        list_parameters.append(function_bond)
                        bondtype_dict[labelf].append(list_parameters)
                    # The label is already present in the bondtype_dict. Checking for any new parameters."
                    else:
                        if labelf in bondtype_dict:
                            lold = labelf
                        if labelr in bondtype_dict:
                            lold = labelr
                        list_parameters = [float(i) for i in tokens[3:]]
                        if len(list_parameters) == 0:
                            is_bond_parameter_duplicates = False
                        else:
                            is_bond_parameter_duplicates = True
                            list_parameters.append(function_bond)
                            bondtype_repeat_dict[lold].append(list_parameters)
                        if is_bond_parameter_duplicates:
                            m = "\n\t\t WARNING: There are paremeters for {} or {} bond\n".format(labelf, labelr)
                            m += "\t\t          in two files: {} and ffbonded.itp\n".format(ff_name)
                            m += "\t\t          Excluding the parameters from {}\n".format(ff_name)
                            m += "\t\t          Bond atoms: {} {} \n".format(iat, jat)
                            m += "\t\t          Please check carefully the parameters.\n".format(ff_name)
                            print(m) if log is None else log.warn(m)
                idx += 1

    # [ angle ] label =============================================================
    is_angles_label = False
    with open(ff_name, 'r') as ff_itp:
        ff_itp.seek(0)
        for num, line in enumerate(ff_itp, 1):
            if re.match(r"^\[(\s*)angles(\s*)]", line):
                d_indices_lines["angles"].append({"file": ff_name, "idx": num - 1})
                is_angles_label = True
        if is_angles_label:
            ff_itp.seek(0)
            lines_itp = ff_itp.readlines()
            idx = d_indices_lines["angles"][0]["idx"] + 1
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
                    function_angle = int(tokens[3])
                    itype = idx_to_atomtype_dict[iat]
                    jtype = idx_to_atomtype_dict[jat]
                    ktype = idx_to_atomtype_dict[kat]
                    labelf = itype + "-" + jtype + "-" + ktype
                    labelr = ktype + "-" + jtype + "-" + itype

                    if labelf not in angletype_dict and \
                            labelr not in angletype_dict:
                        # Function 1,2,5
                        list_parameters = [float(i) for i in tokens[4:]]
                        list_parameters.append(function_bond)
                        angletype_dict[labelf].append(list_parameters)
                    else:
                        if labelf in angletype_dict:
                            lold = labelf
                        if labelr in angletype_dict:
                            lold = labelr
                        list_parameters = [float(i) for i in tokens[4:]]
                        if len(list_parameters) == 0:
                            is_angle_parameter_duplicates = False
                        else:
                            is_angle_parameter_duplicates = True
                            list_parameters.append(function_angle)
                            angletype_repeat_dict[lold].append(list_parameters)
                        if is_angle_parameter_duplicates:
                            m = "\n\t\t WARNING: There are paremeters for {} or {} angle\n".format(labelf, labelr)
                            m += "\t\t          in two files: {} and ffbonded.itp\n".format(ff_name)
                            m += "\t\t          Excluding the parameters from {}\n".format(ff_name)
                            m += "\t\t          Angle atoms: {} {} {} \n".format(iat, jat, kat)
                            m += "\t\t          Please check carefully the parameters.\n".format(ff_name)
                            print(m) if log is None else log.warn(m)
                idx += 1

    # [ dihedral ] label =============================================================
    is_dihedrals_label = False
    with open(ff_name, 'r') as ff_itp:
        ff_itp.seek(0)
        for num, line in enumerate(ff_itp, 1):
            if re.match(r"^\[(\s*)dihedrals(\s*)]", line):
                d_indices_lines["dihedrals"].append({"file": ff_name, "idx": num - 1})
                is_dihedrals_label = True
        if is_dihedrals_label:
            ff_itp.seek(0)
            lines_itp = ff_itp.readlines()
            idx = d_indices_lines["dihedrals"][0]["idx"] + 1
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
                    function_dihedral = int(tokens[4])
                    itype = idx_to_atomtype_dict[iat]
                    jtype = idx_to_atomtype_dict[jat]
                    ktype = idx_to_atomtype_dict[kat]
                    ltype = idx_to_atomtype_dict[lat]
                    labelf = itype + "-" + jtype + "-" + ktype + "-" + ltype
                    labelr = ltype + "-" + ktype + "-" + jtype + "-" + itype
                    # Proper dihedrals
                    if function_dihedral not in [2, 4]:
                        if labelf not in dihedraltype_dict and \
                                labelr not in dihedraltype_dict:
                            list_parameters = [float(i) for i in tokens[5:]]
                            list_parameters.append(function_dihedral)
                            dihedraltype_dict[labelf].append(list_parameters)
                        else:
                            if labelf in dihedraltype_dict:
                                lold = labelf
                            if labelr in dihedraltype_dict:
                                lold = labelr
                            list_parameters = [float(i) for i in tokens[5:]]
                            if len(list_parameters) == 0:
                                is_dihedral_parameter_duplicates = False
                            else:
                                is_dihedral_parameter_duplicates = True
                                dihedraltype_repeat_dict[lold].append(list_parameters)
                            if is_dihedral_parameter_duplicates:
                                m = "\n\t\t WARNING: There are paremeters for {} or {} dihedral\n".format(labelf,
                                                                                                          labelr)
                                m += "\t\t          in two files: {} and ffbonded.itp\n".format(ff_name)
                                m += "\t\t          Excluding the parameters from {}\n".format(ff_name)
                                m += "\t\t          Dihedral atoms: {} {} {} {}\n".format(iat, jat, kat, lat)
                                m += "\t\t          Please check carefully the parameters.\n".format(ff_name)
                                print(m) if log is None else log.warn(m)
                    # Improper dihedrals
                    if function_dihedral in [2, 4]:
                        if labelf not in impropertype_dict and \
                                labelr not in impropertype_dict:
                            list_parameters = [float(i) for i in tokens[5:]]
                            list_parameters.append(function_dihedral)
                            impropertype_dict[labelf].append(list_parameters)
                            improper_idx_list.append([iat, jat, kat, lat])
                        else:
                            if labelf in impropertype_dict:
                                lold = labelf
                            if labelr in impropertype_dict:
                                lold = labelr
                            list_parameters = [float(i) for i in tokens[5:]]
                            if len(list_parameters) == 0:
                                is_improper_parameter_duplicates = False
                            else:
                                is_improper_parameter_duplicates = True
                                impropertype_repeat_dict[lold].append(list_parameters)
                            improper_idx_list.append([iat, jat, kat, lat])
                            if is_improper_parameter_duplicates:
                                m = "\n\t\t WARNING: There are paremeters for {} or {} improper\n".format(labelf,
                                                                                                          labelr)
                                m += "\t\t          in two files: {} and ffbonded.itp\n".format(ff_name)
                                m += "\t\t          Excluding the parameters from {}\n".format(ff_name)
                                m += "\t\t          Improper atoms: {} {} {} {}\n".format(iat, jat, kat, lat)
                                m += "\t\t          Please check carefully the parameters.\n".format(ff_name)
                                print(m) if log is None else log.warn(m)

                idx += 1

    # [ pairs ] label =============================================================
    is_pairs_label = False
    with open(ff_name, 'r') as ff_itp:
        ff_itp.seek(0)
        for num, line in enumerate(ff_itp, 1):
            if re.match(r"^\[(\s*)pairs(\s*)]", line):
                d_indices_lines["pairs"].append({"file": ff_name, "idx": num - 1})
                is_pairs_label = True
        if is_pairs_label:
            ff_itp.seek(0)
            lines_itp = ff_itp.readlines()
            idx = d_indices_lines["pairs"][0]["idx"] + 1
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
                    function_pair = int(tokens[2])
                    itype = idx_to_atomtype_dict[iat]
                    jtype = idx_to_atomtype_dict[jat]
                    labelf = itype + "-" + jtype
                    labelr = jtype + "-" + itype
                    # The label is not in the pairtype_dict, adding the parameters
                    if labelf not in pairtype_dict and \
                            labelr not in pairtype_dict:
                        list_parameters = [float(i) for i in tokens[3:]]
                        list_parameters.append(function_bond)
                        pairtype_dict[labelf].append(list_parameters)
                    # The label is already present in the bondtype_dict. Checking for any new parameters."
                    else:
                        if labelf in pairtype_dict:
                            lold = labelf
                        if labelr in pairtype_dict:
                            lold = labelr
                        list_parameters = [float(i) for i in tokens[3:]]
                        if len(list_parameters) == 0:
                            is_pair_parameter_duplicates = False
                        else:
                            is_pair_parameter_duplicates = True
                            list_parameters.append(function_pair)
                            pairtype_repeat_dict[lold].append(list_parameters)
                        if is_pair_parameter_duplicates:
                            m = "\n\t\t WARNING: There are paremeters for {} or {} pairs\n".format(labelf, labelr)
                            m += "\t\t          in two files: {} and ffbonded.itp\n".format(ff_name)
                            m += "\t\t          Excluding the parameters from {}\n".format(ff_name)
                            m += "\t\t          Bond atoms: {} {} \n".format(iat, jat)
                            m += "\t\t          Please check carefully the parameters.\n".format(ff_name)
                            print(m) if log is None else log.warn(m)
                idx += 1

    return improper_idx_list, bondtype_repeat_dict, angletype_repeat_dict, dihedraltype_repeat_dict, \
        impropertype_repeat_dict, pairtype_repeat_dict
