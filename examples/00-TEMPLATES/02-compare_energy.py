import sys
import os
import re
from collections import defaultdict

# Working directory to extract energies
CHECKDIR=sys.argv[1]
os.chdir(CHECKDIR)
WK=sys.argv[2]
TITLE=sys.argv[3]

dict_gromacs = defaultdict()
dict_lammps = defaultdict()


if  not os.path.isfile("energy.xvg") or not os.path.isfile("clean.log"):
    print("ERROR: Either energy.xvg or clean.log file does not exist")
    exit(2)

# Read energy.xvg from GROMACS
with open("energy.xvg", 'r') as fxvg:
    lines = fxvg.readlines()
    for i in range(0, len(lines)):
        iline = lines[i]
        if re.search(r"^    0.000000  ", iline) is not None:
                gromacs_energyline = iline.split()

with open("energy.xvg", 'r') as fxvg:
    lines = fxvg.readlines()
    for i in range(0, len(lines)):
        iline = lines[i]
        if iline.find("Bond") != -1:
            pos = int(iline.split()[1][1:])+1
            dict_gromacs["Bond"] = float(gromacs_energyline[pos])
        if iline.find("Angle") != -1:
            pos = int(iline.split()[1][1:])+1
            dict_gromacs["Angle"] = float(gromacs_energyline[pos])
        if iline.find("Ryckaert-Bell.") != -1 or iline.find("Proper Dih.") != -1:
            pos = int(iline.split()[1][1:])+1
            try:
                dict_gromacs["Dihedral"] += float(gromacs_energyline[pos])
            except KeyError:
                dict_gromacs["Dihedral"] = float(gromacs_energyline[pos])
        if iline.find("LJ-14") != -1 or iline.find("LJ (SR)") != -1 or iline.find("Disper. corr.") != -1:
            pos = int(iline.split()[1][1:])+1
            try:
                dict_gromacs["Vdw"] += float(gromacs_energyline[pos])
            except KeyError:
                dict_gromacs["Vdw"] = float(gromacs_energyline[pos])

        if iline.find("Coulomb-14") != -1 or iline.find("Coulomb (SR)") != -1 or iline.find("Coul. recip.") != -1:
            pos = int(iline.split()[1][1:])+1
            try:
                dict_gromacs["Coulomb"] += float(gromacs_energyline[pos])
            except KeyError:
                dict_gromacs["Coulomb"] = float(gromacs_energyline[pos])
        if iline.find("Improper") != -1:
            pos = int(iline.split()[1][1:])+1
            dict_gromacs["Improper"] = float(gromacs_energyline[pos])

# Read clean.log from LAMMPS
with open("clean.log", 'r') as flog:
    lines = flog.readlines()
    for i in range(0, len(lines)):
        iline = lines[i]
        if re.search(r"^         0   0", iline) is not None:
            lammps_energyline = iline.split()

with open("clean.log", 'r') as flog:
    lines = flog.readlines()
    for i in range(0, len(lines)):
        iline = lines[i]
        if iline.find("E_bond") != -1:
            tokens = iline.split()
            dict_lammps["Bond"] = float(lammps_energyline[int(tokens.index("E_bond"))])
            dict_lammps["Angle"] = float(lammps_energyline[int(tokens.index("E_angle"))])
            dict_lammps["Dihedral"] = float(lammps_energyline[int(tokens.index("E_dihed"))])
            dict_lammps["Vdw"] = float(lammps_energyline[int(tokens.index("E_vdwl"))])
            dict_lammps["Coulomb"] = float(lammps_energyline[int(tokens.index("E_coul"))])
            dict_lammps["Coulomb"] += float(lammps_energyline[int(tokens.index("E_long"))])
            dict_lammps["Improper"] = float(lammps_energyline[int(tokens.index("E_impro"))])

print (dict_gromacs)
print(dict_lammps)


# Write results
logname = os.path.join(WK, "out_tests.log")
keys = ["Bond", "Angle", "Dihedral", "Vdw", "Coulomb", "Improper"]
with open(logname, 'a') as flog:
    line = "\n{0:>30s} {1:>18s} {2:>14s} {3:>10s}\n".format("GROMACS", "LAMMPS", "Diff", "(in kcal/mol)")
    flog.writelines(line)
    for key in keys:

        try:
            e1 = float(dict_gromacs[key])/4.184
            e2 = float(dict_lammps[key])
            line = "{0:>14s} {1:18.7f} {2:18.7f} {3:18.7f} ({4:10.5f}%)\n".format(key, e1, e2, e1-e2, e1/e2*100.)
        except ZeroDivisionError:
            e1 = float(dict_gromacs[key])/4.184
            e2 = float(dict_lammps[key])
            line = "{0:>14s} {1:18.7f} {2:18.7f} {3:18.7f} ({4:10.5f}%)\n".format(key, e1, e2, e1-e2, 100.)
        except KeyError:
            line = "{0:>14s} {1:18.7f} {2:18.7f} {3:18.7f} ({4:10.5f}%)\n".format(key, 0.0, 0.0, 0.0, 100.)
        flog.writelines(line)






exit(0) 