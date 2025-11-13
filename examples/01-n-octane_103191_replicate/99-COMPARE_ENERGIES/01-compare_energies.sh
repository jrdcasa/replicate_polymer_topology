#!/bin/bash

extract_energy_opls () {

    # GROMACS ============================================================
    GROENE00="../99-GROMACS_INP_MIN_NVT/energy.xvg"
    GROLINE=`egrep "    0\.0000" $GROENE00`
    EBONDGRO=`echo $GROLINE | awk '{printf "%.4f", $2}'`    #kJ/mol
    EANGLGRO=`echo $GROLINE | awk '{printf "%.4f", $3}'`    #kJ/mol
    ETORSGRO=`echo $GROLINE | awk '{printf "%.4f", $4}'`    #kJ/mol
    EVDWGRO=`echo $GROLINE  | awk '{printf "%.4f", $7+$5}'`    #kJ/mol
    ECOULGRO=`echo $GROLINE | awk '{printf "%.4f", $6+$8+$9}'`    #kJ/mol
    EPOTGRO=`echo $GROLINE  | awk '{printf "%.4f", $10}'`    #kJ/mol

    # LAMMPS =============================================================
    LMPENE00="../99-LAMMPS_INP_MIN_NVT/log.lammps"
    LMPLINE=`egrep "^         0" $LMPENE00`
    EBONDLMP=`echo $LMPLINE | awk '{print $3}'`    #kcal/mol
    EANGLLMP=`echo $LMPLINE | awk '{print $4}'`    #kcal/mol
    ETORSLMP=`echo $LMPLINE | awk '{print $5}'`    #kcal/mol
    EVDWLMP=`echo $LMPLINE  | awk '{printf "%.4f", $8}'`    #kJ/mol
    ECOULLMP=`echo $LMPLINE | awk '{print $9+$10}'`    #kJ/mol
    EPOTLMP=`echo $LMPLINE  | awk '{printf "%.4f", $12}'`    #kJ/mol

    # CALCULATIONS =======================================================
    DELTAEBOND=`echo "$EBONDGRO $EBONDLMP"|awk '{printf "%.4f", $1-$2*4.184}'`
    DELTAEANGL=`echo "$EANGLGRO $EANGLLMP"|awk '{printf "%.4f", $1-$2*4.184}'`
    DELTAETORS=`echo "$ETORSGRO $ETORSLMP"|awk '{printf "%.4f", $1-$2*4.184}'`
    DELTAEIMP=`echo "$EIMPGRO $EIMPLMP"|awk '{printf "%.4f", $1-$2*4.184}'`
    DELTAECOUL=`echo "$ECOULGRO $ECOULLMP"|awk '{printf "%.4f", $1-$2*4.184}'`
    DELTAEVDW=`echo "$EVDWGRO $EVDWLMP" |awk '{printf "%.4f", $1-$2*4.184}'`
    DELTAEPOT=`echo "$EPOTGRO $EPOTLMP" |awk '{printf "%.4f", $1-$2*4.184}'`

    echo "                                                  GROMACS-LAMMPS"
    echo "Ebond(GROMACS, kJ/mol)   Ebond(LAMMPS, kcal/mol)   Edelta(kJ/mol)"
    echo "    $EBONDGRO                    $EBONDLMP           $DELTAEBOND"
    echo ""
    echo "Ebend(GROMACS, kJ/mol)   Ebend(LAMMPS, kcal/mol)   Edelta(kJ/mol)"
    echo "    $EANGLGRO                    $EANGLLMP           $DELTAEANGL"
    echo ""
    echo "Etors(GROMACS, kJ/mol)   Etors(LAMMPS, kcal/mol)   Edelta(kJ/mol)"
    echo "    $ETORSGRO                    $ETORSLMP           $DELTAETORS"
    echo ""
#    echo "Eimp(GROMACS, kJ/mol)   Eimp(LAMMPS, kcal/mol)   Edelta(kJ/mol)"
#    echo "    $EIMPGRO                    $EIMPLMP           $DELTAEIMP"
#    echo ""
    echo "ECoul(GROMACS, kJ/mol)   ECoul(LAMMPS, kcal/mol)   Edelta(kJ/mol)"
    echo "    $ECOULGRO                    $ECOULLMP           $DELTAECOUL"
    echo ""
    echo "EVdW(GROMACS, kJ/mol)   EVdW(LAMMPS, kcal/mol)   Edelta(kJ/mol)"
    echo "    $EVDWGRO                 $EVDWLMP         $DELTAEVDW"
    echo ""
    echo "EPOT(GROMACS, kJ/mol)   EPOT(LAMMPS, kcal/mol)   Edelta(kJ/mol)"
    echo "    $EPOTGRO                 $EPOTLMP         $DELTAEPOT"
    
    echo $GROLINE | awk '{print $5+$6+$7+$8+$9+$10+$2+$3+$4}'
    echo $LMPLINE | awk '{print $3+$4+$5+$6+$9+$8+$10}'

}



# =========================== MAIN ==========================
WK=`pwd`
# Example 01
extract_energy_opls
#cd ${WK}
