# REPLICATE_POLYMER: A tool to replicate periodically a molecule or polymer chain

`replicate_polymer` allows one to create a molecular system by replicating a molecule in a simulation box along with the application of a force field to the atoms. The program produces input files to be used in GROMACS and LAMMPS.  

## Basic usage

After installation you can simply write **replicate_polymer** to run the program.

```bash
$ replicate_polymer -h
Warning: importing 'simtk.openmm' is deprecated.  Import 'openmm' instead.
usage: replicate_polymer.py [-h] -p PDB_FILE -f XML_FILE --images image_x image_y image_z [-e MDENGINE] [--noh | --index INDEX] [--boxlength a b c]
                            [--boxangle alpha beta gamma]

Replicate and apply a force field to a polymer or molecule.

optional arguments:
  -h, --help            show this help message and exit
  -p PDB_FILE, --pdb PDB_FILE
                        A pdb file containing the structure to be replicated.
  -f XML_FILE, --forcefield XML_FILE
                        A XML file containing the force field parameters using the FOYER format
  --images image_x image_y image_z
                        A list with the number of images to replicate in the x, y and z axis.
  -e MDENGINE, --engine MDENGINE
                        MD package to perform the calculations.
  --noh                 Remove hydrogens for a united atom representation.
  --index INDEX         Indices of atoms to be removed from the PDB.
  --boxlength a b c     Box lengths in nanometers.
  --boxangle alpha beta gamma
                        Box angles in degrees.
```

### [Installation instructions](doc/installation.md)

The best way to learn how this program works is through examples:

### [Example 1: Creation of a n-octane crystal](doc/example01.md)

### [Example 2: Creation of a simulation box with a polymer chain](doc/example02.md)

### [Example 3: Creation of a folded crystal](doc/example03.md)
