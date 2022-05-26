
 -- C150_1fold_initial.pdb ----------------------------------------------------
 Preparado a mano en MS usando el cristal CIF del PE. Primero se crea la supercelda para que tenga 150 C y despues se construye una estructura no-peridoca. Finalmente se exporta a pdb

 -- C150_1fold_initial_new.pdb ----------------------------------------------------

Ordenar atomos. He usado gaussview.
CRYST1    7.388    4.929   91.404  90.00  90.00  90.00 PANAM          1

 -- 10x10x1 Replicate ----------------------------------------------------
home/jramos/PycharmProjects/MOSDEF_SOFTWARE/MyPolymers/topocordpolymer.py
  -p C150_1fold_initial_new.pdb -f ../oplsaa.xml --images 10 10 1 -e lammps
