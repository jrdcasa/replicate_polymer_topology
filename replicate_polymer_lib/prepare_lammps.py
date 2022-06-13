import os
import datetime
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
    ll = LammpsParser(lammps_template_input, system=system)
    ll.write(nonbonded_style=nonbonded_style_defaults)
