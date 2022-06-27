#!bin/bash

prepare_simulations() {

   # Copy gromacs files
   cp  $1/*_replicate.gro $1/*_replicate.top $2
   # Copy lammps intermol files
   cp  $1/*_replicate.lmp $1/*_replicate.inp $2
   # Copy lammps clean files
   cp  $1/*_replicate_clean.lmp $1/*_replicate_clean.inp $2
   # Copy mdp
   cp $1/minim_template.mdp $2

   cd $2
   # RUN GROMACS
   module load GROMACS
   gmx grompp -f minim_template.mdp -c `ls *_replicate.gro` -p `ls *.top` -o new.tpr -maxwarn 10
   gmx mdrun -s new.tpr
   echo 1 2 3 4 5 6 7 8 9 10 0 | gmx energy -f ener.edr

   rm \#*
   # RUN LAMMPS
   module purge
   module load lammps

   # Clean
   INP=`ls *_clean.inp`
   lmp_mpi -in $INP -log clean.log
   
   # Intermol
   echo "===================="
   echo "Modify and RUN LAMMPS $2: `ls ${2}/*_replicate.inp` `ls ${2}/*_replicate.lmp`"
   echo ""
   echo "cd  `ls ${2}/*_replicate.inp`"
   echo ""
   echo "lmp_mpi -in  `ls ${2}/*_replicate.inp` -log intermol.log"
   echo "===================="

}



# ========================== MAIN =====================================
WK=`pwd`

source /home/jramos/Programacion/sandboxes/sandbox_common/bin/activate
mkdir -p XX-CHECK_ENERGY

# Example 01
TESTDIR=01-n-octane_103191_replicate
cd ${TESTDIR}
replicate_polymer -p 103191_noctane_order_cryst.pdb -f ../../forcefields/oplsaa.xml --images 6 6 6 --engine lammps
mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
cp ${WK}/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
cd $WK
