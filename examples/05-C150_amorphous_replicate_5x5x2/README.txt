mkdir 01-CREATE
cd 01-CREATE
# Modify the simulation box CRYST1
replicate_polymer -p C150_1fold_initial_new.pdb -f ../../../forcefields/oplsaa.xml --images 1 1 1 --engine lammps

cd ..
mkdir 02-LAMMPS_amorphous
cd 02-LAMMPS_amorphous
# Modify inp file
lmp_pysimm -in C150_1fold_initial_new_replicate.inp

cd ..
mkdir 03-REPLICATE_AMORPHOUS
# Extract the last structure to pdb with connect table

