#!/bin/bash
#SBATCH --partition=gpu_rtx3090_14p
#SBATCH -n 14
#SBATCH -N 1
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=3G
#SBATCH --job-name=DVE_OPT

#module purge
#module load foss GROMACS/2021.5
source /usr/bin/GMXRC.bash


MDP="eq2_NPT.mdp"
TPR="../02-MD_EQ01_NVT/new_topo.tpr"
CPI="../02-MD_EQ01_NVT/state.cpt"
TOP="../00-FORCE_FIELD_POLYPLY/PVL_OPLS.top"

gmx grompp -f $MDP -p $TOP -c $TPR -t $CPI -o new_topo.tpr --maxwarn 5 >& outgro.dat
gmx mdrun -gpu_id 0 -pin auto -nt 14 -bonded gpu -update gpu -nb gpu -s new_topo.tpr -v -noappend  >& mdout.dat
#gmx mdrun -s new_topo.tpr -v -noappend  >& outmd.dat
#echo 0|gmx trjconv -f traj_comp.part0001.xtc -o traj_unwrap.xtc -pbc whole -s new_topo.tpr
rm \#*\n

