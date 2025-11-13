#!bin/bash

prepare_simulations() {

   # Copy gromacs files
   cp  $1/*_replicate.gro $1/*_replicate.top $2
   # Copy lammps clean files
   cp  $1/*_replicate_clean.lmp $1/*_replicate_clean.inp $2
   # Copy mdp
   #cp $1/minim_template.mdp $2

   cd $2
   # RUN GROMACS
   module load GROMACS
   gmx grompp -f minim_template.mdp -c `ls *_replicate.gro` -p `ls *.top` -o new.tpr -maxwarn 10
   gmx mdrun -s new.tpr
   echo 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 0 | gmx energy -f ener.edr

   rm \#*
   # RUN LAMMPS
   module purge
   module load lammps
   INP=`ls *_clean.inp`
   lmp_mpi -in $INP -log clean.log
   
}

compare_energy() {

   python ${TEMPLATES}/02-compare_energy.py $1 $2 $3
   EXITCODE=$?
   if [[ ${EXITCODE} -ne 0 ]]; then
      echo "Example $3: Energies cannot be compare either energy.xvg or clean.log are not correct." >>$2/out_tests.log
   fi

}

# ========================== MAIN =====================================
WK=`pwd`
rm -r XX-CHECK_ENERGY
cd 00-TEMPLATES
TEMPLATES=`pwd`
cd $WK
if [[ $# -eq 0 ]]; then
    runjob=0        # Run all examples 
else
    runjob=$1  
fi

module purge

echo "========== TEST EXAMPLES ==========" >${WK}/out_tests.log

source /home/jramos/Programacion/sandboxes/sandbox_common/bin/activate
mkdir -p XX-CHECK_ENERGY

# Example 01 ************************************
if [[ $runjob -eq 0 || $runjob -eq 1 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 01: 01-n-octane_6x6x6" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=01-n-octane_103191_replicate
    cd ${TESTDIR}
    replicate_polymer -p 103191_noctane_order_cryst.pdb -f ../../forcefields/oplsaa.xml --images 6 5 5 -e lammps
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 1
    cd $WK
fi
# Example 01 ************************************

# Example 02 ************************************
if [[ $runjob -eq 0 || $runjob -eq 2 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 02: 02-C150_amorphous_singlechain" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=02-C150_amorphous_singlechain
    cd ${TESTDIR}
    replicate_polymer -p C150_amorphous_residues_boxinfo.pdb -f ../../forcefields/oplsaa.xml --images 5 5 5 --engine lammps --boxlength 8.0 8.0 3.0 --boxangle 65.0 60.0 70.0   
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 2
    cd $WK
fi
# Example 02 ************************************

# Example 03 ************************************
if [[ $runjob -eq 0 || $runjob -eq 3 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 03: 03-C150_onefolded_replicate_5x5x5" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=03-C150_onefolded_replicate_5x5x5
    cd ${TESTDIR}
    replicate_polymer -p C150_1fold_residues.pdb -f ../../forcefields/oplsaa.xml --images 5 5 5 -e lammps   
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 3
    cd $WK
fi
# Example 03 ************************************

# Example 04 ************************************
if [[ $runjob -eq 0 || $runjob -eq 4 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 04: 04-Trimer_P4HB_OPLSAA" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=04-Trimer_P4HB_OPLSAA
    cd ${TESTDIR}
    replicate_polymer -p TrimerP4HB_residues.pdb -f ../../forcefields/oplsaa.xml --images 2 2 2 -e lammps --boxlength 2.3 2.3 2.3    
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 4
    cd $WK
fi
# Example 04 ************************************

# Example 05 ************************************
if [[ $runjob -eq 0 || $runjob -eq 5 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 05: 05-Trimer_P4HB_TrappeUA" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=05-Trimer_P4HB_TrappeUA
    cd ${TESTDIR}
	replicate_polymer -p TrimerP4HB_residues.pdb -f ../../forcefields/trappe-ua.xml --images 1 1 1 -e lammps --boxlength 2.3 2.3 2.3 --noh --verbose
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 5
    cd $WK
fi
# Example 05 ************************************

# Example 06 ************************************
if [[ $runjob -eq 0 || $runjob -eq 6 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 06: 06-pdb_two_diferent_molecules" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=06-pdb_two_diferent_molecules
    cd ${TESTDIR}
	replicate_polymer -p two_molecules.pdb -f ../../forcefields/oplsaa.xml --images 2 2 2 -e lammps --boxlength 2.3 2.3 2.3
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 6
    cd $WK
fi
# Example 06 ************************************

# Example 07 ************************************
if [[ $runjob -eq 0 || $runjob -eq 7 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 07: 07-blend_xsd" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=07-blend_xsd
    cd ${TESTDIR}
    replicate_polymer -p blend_residues.pdb -f ../../forcefields/oplsaa.xml --images 2 2 2 -e lammps
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 7
    cd $WK
fi
# Example 07 ************************************

# Example 08 ************************************
if [[ $runjob -eq 0 || $runjob -eq 8 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 08: MIX PE + PE-SCB OPLS" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=08-blend_toxwaerd1mol_OPLS
    cd ${TESTDIR}
    replicate_polymer -p ./01-PREPARE/Toxwaerd_twokinds_residues.pdb -f ../../forcefields/oplsaa.xml --images 2 2 3 -e lammps --boxlength 2.3 2.3 2.3
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 8
    cd $WK
fi
# Example 08 ************************************

# Example 09 ************************************
if [[ $runjob -eq 0 || $runjob -eq 9 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 09: MIX PE + PE-SCB TrappeUA-Toxwaerd (No improper)" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=09-blend_toxward1mol_TrappeTox
    cd ${TESTDIR}
    replicate_polymer -p ./01-PREPARE/Toxwaerd_twokinds_residues.pdb -f ../../forcefields/trappe-ua_PEToxvaerd.xml --images 2 2 3 -e lammps --noh --boxlength 2.3 2.3 2.3
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 9
    cd $WK
fi
# Example 09 ************************************

# Example 10 ************************************
if [[ $runjob -eq 0 || $runjob -eq 10 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 10: MIX PE + PE-SCB TrappeUA-Toxwaerd (Improper file)" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=10-blend_toxward1mol_TrappeTox_improperfile
    cd ${TESTDIR}
    replicate_polymer -p ./01-PREPARE/Toxwaerd_twokinds_residues.pdb -f ../../forcefields/trappe-ua_PEToxvaerd.xml --images 2 2 3 -e lammps --noh --boxlength 2.3 2.3 2.3 --impropers improper.ndx
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 10
    cd $WK
fi
# Example 10 ************************************

# Example 11 ************************************
if [[ $runjob -eq 0 || $runjob -eq 11 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 11: MIX PE + PE-SCB LOPLS No consequtive" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=11-blend_noconsecutive
    cd ${TESTDIR}
    replicate_polymer -p blend2.pdb -f ../../forcefields/loplsaa.xml -e lammps --images 2 2 2
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 11
    cd $WK
fi
# Example 11 ************************************
# Example 12 ************************************
if [[ $runjob -eq 0 || $runjob -eq 12 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 12: Simulation ETH-EVA 50%-wt random copolymer" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=12-ETH-EVA_50wt_trappeua
    cd ${TESTDIR}
    replicate_polymer -p ETH-EVA_50wt_0000.pdb -f ../../forcefields/trappe-ua.xml -e lammps --images 1 1 1 --noh --verbose
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 12
    cd $WK
fi
# Example 12 ************************************
# Example 13 ************************************
if [[ $runjob -eq 0 || $runjob -eq 13 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 13: Simulation ETH-EVA 50%-wt random copolymer" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=13-ETH-EVA_50wt_oplsaa
    cd ${TESTDIR}
    replicate_polymer -p ETH-EVA_50wt_0000.pdb -f ../../forcefields/oplsaa.xml -e lammps --images 1 1 1 
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 13
    cd $WK
fi
# Example 13 ************************************

# Example 23 ************************************
if [[ $runjob -eq 0 || $runjob -eq 23 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 23: Single chain isotactic PP (OPLS)" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=23-iPP_singleChain_OPLS
    cd ${TESTDIR}
    replicate_polymer -p iPP_SC_40mon_50Ch_residues_center.pdb -f ../../forcefields/oplsaa.xml --images 2 2 2 -e lammps
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 23
    cd $WK
fi
# Example 23 ************************************
# Example 24 ************************************
if [[ $runjob -eq 0 || $runjob -eq 24 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 24: Single chain isotactic PP (Trappe-UA) improper guess" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=24-iPP_singleChain_TrappeUA_imprguess
    cd ${TESTDIR}
    replicate_polymer -p iPP_SC_40mon_50Ch_residues_center.pdb -f ../../forcefields/trappe-ua.xml --images 2 2 2 -e lammps --noh --impropers guess
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 24
    cd $WK
fi
# Example 24 ************************************
# Example 25 ************************************
if [[ $runjob -eq 0 || $runjob -eq 25 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 25: 50Chains 40Mon isotactic PP (OPLS)" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=25-iPP_AC_OPLS
    cd ${TESTDIR}
    replicate_polymer -p iPP_SC_40mon_50Ch_residues.pdb -f ../../forcefields/oplsaa.xml --images 1 1 1 -e lammps
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 25
    cd $WK
fi
# Example 25 ************************************
# Example 26 ************************************
if [[ $runjob -eq 0 || $runjob -eq 26 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 26: 50Chains 40Mon isotactic PP (TRAPE-UA) improper file" >>${WK}/out_tests.log
    echo "Impropers are not all iso. This is for the MS initial strucuture"  >>${WK}/out_tests.log
    echo "I have to investigate why????"  >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=26-iPP_AC_TrappeUA_imprfile
    cd ${TESTDIR}
    replicate_polymer -p iPP_SC_40mon_50Ch_residues.pdb -f ../../forcefields/trappe-ua.xml --images 1 1 1 -e lammps  --noh --impropers improper.ndx
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 26
    cd $WK
fi
# Example 26 ************************************
# Example 27 ************************************
if [[ $runjob -eq 0 || $runjob -eq 27 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 27: Single Chain isotactic PP (TRAPE-UA) improper file" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=27-iPP_singleChain_TrappeUA_imprfile
    cd ${TESTDIR}
    replicate_polymer -p iPP_SC_40mon_50Ch_residues_center.pdb -f ../../forcefields/trappe-ua.xml --images 2 2 2 -e lammps --noh --impropers improper.ndx
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 27
    cd $WK
fi
# Example 27 ************************************
# Example 28 ************************************
if [[ $runjob -eq 0 || $runjob -eq 28 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 28: Test" >>${WK}/out_tests.log
    echo "Impropers are not all iso. This is for the MS initial structure"  >>${WK}/out_tests.log
    echo "I have to investigate why????"  >>${WK}/out_tests.log
    echo "In this case the dihedral 50       49       52       51  2 is positive instead of negative" >>${WK}/out_tests.log
    echo "I have to investigate why????"  >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=28-TEST
    cd ${TESTDIR}
    replicate_polymer -p 0123.pdb -f ../../forcefields/trappe-ua.xml --images 1 1 1 -e lammps --noh --impropers improper.ndx
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 28
    cd $WK
fi 
# Example 28 ************************************
# Example 29 ************************************
if [[ $runjob -eq 0 || $runjob -eq 29 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 29: Blend iPP-PE 25-25Ch 40monomers each polymer OPLS" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=29-PE-iPP_AC_50Ch_OPLS
    cd ${TESTDIR}
    replicate_polymer -p iPP-PE_40mon_25-25Ch_residues.pdb -f ../../forcefields/oplsaa.xml --images 1 1 1 -e lammps
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 29
    cd $WK
fi
# Example 29 ************************************
# Example 30 ************************************
if [[ $runjob -eq 0 || $runjob -eq 30 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 30: Blend iPP-PE 25-25Ch 40monomers each polymer LOPLS" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=30-PE-iPP_AC_50Ch_LOPLS
    cd ${TESTDIR}
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 30
    cd $WK
fi
# Example 30 ************************************
# Example 31 ************************************
if [[ $runjob -eq 0 || $runjob -eq 31 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 31: Blend iPP-PE 25-25Ch 40monomers each polymer TrappeUA" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=31-PE-iPP_AC_50Ch_Trappe
    cd ${TESTDIR}
    replicate_polymer -p iPP-PE_40mon_25-25Ch_residues.pdb -f ../../forcefields/trappe-ua.xml --images 1 1 1 -e lammps --noh --impropers improper.ndx
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 31
    cd $WK
fi
# Example 31 ************************************
# Example 32 ************************************
if [[ $runjob -eq 0 || $runjob -eq 32 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 32: Blend iPP-PE 25-25Ch 40monomers each polymer TrappeUATox" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=32-PE-iPP_AC_50Ch_TrappeTox
    cd ${TESTDIR}
    replicate_polymer -p iPP-PE_40mon_25-25Ch_residues.pdb -f ../../forcefields/trappe-ua_PEToxvaerd.xml --images 1 1 1 -e lammps --noh --impropers improper.ndx
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 32
    cd $WK
fi
# Example 32 ************************************
# Example 33 ************************************
if [[ $runjob -eq 0 || $runjob -eq 33 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 33: Blend iPP-PE 25-25Ch 40monomers (2x2x2) each polymer LOPLS/OPLS" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=33-PE-iPP_AC_50Ch_LOPLS_2_2_2
    cd ${TESTDIR}
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 33
    cd $WK
fi
# Example 33 ************************************
# Example 34 ************************************
if [[ $runjob -eq 0 || $runjob -eq 34 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 34: Blend iPP-PE 25-25Ch 40monomers (2x2x2) each polymer LOPLS/OPLS" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=34-PE-iPP_AC_50Ch_TrappeTox_2_2_2
    cd ${TESTDIR}
    replicate_polymer -p iPP-PE_40mon_25-25Ch_residues.pdb -f ../../forcefields/trappe-ua_PEToxvaerd.xml --images 2 2 2 -e lammps --noh --impropers improper.ndx
    mkdir -p ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    cp ${WK}/00-TEMPLATES/minim_template.mdp ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    prepare_simulations ${WK}/${TESTDIR} ${WK}/XX-CHECK_ENERGY/${TESTDIR}
    compare_energy ${WK}/XX-CHECK_ENERGY/${TESTDIR} ${WK} 34
    cd $WK
fi
# Example 33 ************************************
