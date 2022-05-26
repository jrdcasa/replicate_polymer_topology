import os
from intermol.convert import _load_gromacs
from intermol.lammps.lammps_parser import LammpsParser


# =============================================================================
def prepare_lammps(basefile, title=None, nonbonded_style_defaults=None):

    if nonbonded_style_defaults is None:
        nonbonded_style_defaults = 'pair_style lj/cut/coul/long 9.0 9.0\nkspace_style pppm 1e-6\n'

    system, prefix, gro_in, top_in = _load_gromacs(basefile)
    if title is not None:
        system.name = title
    lammps_template_input = os.path.splitext(basefile[0])[0] + ".inp"
    ll = LammpsParser(lammps_template_input, system=system)
    ll.write(nonbonded_style=nonbonded_style_defaults)
