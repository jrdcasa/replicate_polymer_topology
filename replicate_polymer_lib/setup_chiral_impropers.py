from _collections import defaultdict
from ext_libc.c_distC_replicate import dihedralC
import re


# =============================================================================
def setup_chiral_impropers(filenamepdb, improper_filename=None, logging=None):

    """
    Get improper angles from a file. If filename is None then the improper are guessed
    but this only works with alkanes. Otherwise, the number of impropers returns is zero

    Args:
        filenamepdb:
        improper_filename:
        logging:

    Returns:
        A tuple [number of impropers, list of impropers in GROMACS format]

    """

    # Create a dictionary with CONNECT section. Classify backbone atoms
    graph = defaultdict()
    is_backbone = []
    xyz = []
    element = []
    chain = []
    nimpropers_list = 0
    improper_list = []
    nimpropers = 0
    with open(filenamepdb, 'r') as fpdb:
        lines = fpdb.readlines()
        for iline in lines:
            if iline.find("ATOM") != -1 or iline.find("HETATM") != -1:
                try:
                    # backbone = 0.00, not backbone = 1.00
                    is_backbone.append(int(float(iline[60:66])))
                except ValueError:
                    is_backbone.append(int(1))
                xyz.append([float(iline[30:38]),
                            float(iline[38:46]),
                            float(iline[46:54])])
                element.append(iline[76:78])
                chain.append(iline[21:22])
            if iline.find("CONECT") != -1:

                ll = [iline[0:6]]
                idx = 6
                while idx < len(iline):
                    ll.append(iline[idx:idx+5])
                    idx += 5
                # print(iline, len(iline), ll)
                key = int(ll[1])-1
                values = [int(i)-1 for i in ll[2:-1]]
                graph[key] = values

    # Try to guess the improper angles. The rule is as follow:
    #   1. A branch point has 3 or more neigbours and it must be a backbone atom.
    #   2. There are so many impropers as branch atoms in the neighbours.
    #   3. Improper angles on the branches are not taken into account.
    #   4. Atom 1 is the branch point on the backbone chain
    #   5. Atom 4 is the branch atom
    #   6. Atom 2 and 3 are in order
    #    Example: (5-4-6-10)
    #              HEAD -----  (4)--(5)--(6) ----- TAIL
    #                                |
    #                              (10)

    improper_lines = "; improper angles\n"
    improper_lines += ";  ai    aj    ak    al funct  c0  c1  c2  c3  c4  c5\n"
    if improper_filename is None:
        # nkinds = len(set(chain))
        # improper_list = []
        # nimpropers_list = []
        for key, values in graph.items():
            if len(values) < 3:
                continue
            # if is_backbone[key] == 0:
            if is_backbone[key] == 1:
                continue
            atoms_improper = [key]
            br_onwait = []
            v = sorted(values)
            for iat in v:
                if is_backbone[iat] == 1:
                    atoms_improper.append(iat)
                else:
                    if element[iat].count('C') != 0:
                        br_onwait.append(iat)
            for iat in br_onwait:
                atoms_improper.append(iat)

            try:
                at1, at2, at3, at4 = atoms_improper[0:4]
                d = dihedralC(xyz[at1][0], xyz[at1][1], xyz[at1][2],
                              xyz[at2][0], xyz[at2][1], xyz[at2][2],
                              xyz[at3][0], xyz[at3][1], xyz[at3][2],
                              xyz[at4][0], xyz[at4][1], xyz[at4][2],)
                nimpropers += 1
                if d < 0:
                    dd = -30.5
                else:
                    dd = 30.5

                # ibbl = []
                # ibrl = []
                # ibpl = [atoms_improper[0]]
                # for iat in atoms_improper[1:]:
                #     if is_backbone[iat] == 0:
                #         ibbl.append(iat)
                #     else:
                #         ibrl.append(iat)
                # ibbl = sorted(ibbl)
                # at1 = ibpl[0]
                # at2 = ibbl[0]
                # at3 = ibbl[1]
                # at4 = ibrl[0]
                improper_lines += "{0:8d} {1:8d} {2:8d} {3:8d}  2   {4:7.3f} 0.5178E+03\n".\
                    format(at1+1, at2+1, at3+1, at4+1, dd)
            except ValueError:
                pass
            nimpropers_list = [nimpropers]
            improper_list = [improper_lines]
    else:
        with open(improper_filename, 'r') as fimp:
            lines = fimp.readlines()
            # nkinds = int(lines[0].split()[1])
            improper_list = []
            nimpropers_list = []
            idx = 1
            while idx < len(lines):
                if re.match("^#", lines[idx]):
                    idx += 1
                    continue

                a, b = lines[idx][:-1].split()
                b = int(b)
                nimpropers_list.append(b)
                idx += 1
                improper_lines = ''
                for j in range(idx, idx+b):
                    iline = lines[j].split()
                    if int(iline[4]) == 4:
                        improper_lines += "{0:8d} {1:8d} {2:8d} {3:8d}  {4:1d}   {5:7.3f} {6:s} {7:d}\n".\
                                          format(int(iline[0]), int(iline[1]), int(iline[2]), int(iline[3]),
                                                 int(iline[4]), float(iline[5]), iline[6], int(iline[7]))
                    elif int(iline[4]) == 2:
                        improper_lines += "{0:8d} {1:8d} {2:8d} {3:8d}  {4:1d}   {5:7.3f} {6:s}\n".\
                                          format(int(iline[0]), int(iline[1]), int(iline[2]), int(iline[3]),
                                                 int(iline[4]), float(iline[5]), iline[6])
                    else:
                        m = "\n\t\tERROR: Improper type {} is not implemented".format(int(iline[4]))
                        logging.error(m)
                improper_list.append(improper_lines)
                idx += b

    return nimpropers_list, improper_list


# =============================================================================
def improper_test():

    #         O(l)
    #         ||
    # O(k) - C(j) - C(i)
    ci = [1.749, -5.485, 1.191]
    cj = [0.981, -4.161, 1.349]
    ck = [-0.486, -4.175, 1.709]
    cl = [1.695, -2.846, 1.145]

    d1 = dihedralC(ci[0], ci[1], ci[2],
                   cj[0], cj[1], cj[2],
                   ck[0], ck[1], ck[2],
                   cl[0], cl[1], cl[2],)

    print(d1)
