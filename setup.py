import os
import glob
from setuptools import setup, find_packages, Extension
import numpy

# -----------------------------------------------------------------------------
# Read version
# -----------------------------------------------------------------------------
version = {}
try:
    with open("version.py") as f:
        exec(f.read(), version)
except FileNotFoundError:
    version["__version__"] = "N/A"


# Try to import Cython — fall back to .c files if not available
try:
    from Cython.Build import cythonize
    use_cython = True
except ImportError:
    def cythonize(*args, **kwargs):
        pass
    use_cython = False

# Find .pyx files inside topology/ext_libc
extlib_dir = "replicate_polymer_lib/ext_libc"
pyx_files = glob.glob(os.path.join(extlib_dir, "*.pyx"))

# If Cython is not installed, use the corresponding .c files
if use_cython:
    sources = pyx_files
else:
    sources = [f.replace(".pyx", ".c") for f in pyx_files]

# Build one Extension object per .pyx/.c module
extensions = [
    Extension(
        name=f"replicate_polymer_lib.ext_libc.{os.path.splitext(os.path.basename(src))[0]}",
        sources=[src],
        include_dirs=[numpy.get_include()],
        extra_compile_args=["-O3", "-fopenmp"],
        extra_link_args=["-fopenmp"],
    )
    for src in sources
]

# Cythonize if available
if use_cython:
    ext_modules = cythonize(
        extensions,
        compiler_directives={
            "language_level": "3",
            "boundscheck": False,
            "wraparound": False,
            "cdivision": True,
        },
        annotate=False,
)
else:
    ext_modules = extensions

setup(
    name="replicate_polymer",
    version=version["__version__"],
    author="Javier Ramos",
    author_email="jrdcasa@gmail.com",
    description="It allows one to create a molecular system in a simulation box with the application of "
                "a force field to the atoms. The program produces input files to be used in GROMACS and LAMMPS.",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/jrdcasa/topology",
    license="GPL-3.0-or-later",
    packages=find_packages()+["replicate_polymer_lib.ext_libc"],
    include_package_data=True,
    install_requires=[
        "argcomplete",
        "mdtraj",
        "openmm",
        "parmed",
        "lark",
        "ele",
        "lxml"
    ],
    python_requires=">=3.6",
    entry_points={
        "console_scripts": [
            "replicate_polymer = replicate_polymer:main",
            "create_library_monomer = create_library_monomer:main",
            "create_lammps_from_gromacs = create_lammps_from_gromacs:main",
            "top2forcefieldxml = top2forcefieldxml:main",
            "remove_hydrogens = remove_hydrogens:main"
        ]
    },
    ext_modules=ext_modules,
    # cmdclass={
    #     "build_ext": BuildIndigox,
    # },
)
