# Minimization
gmx grompp -f minim.mdp -c ../103191_noctane_order_cryst_replicate.gro -p ../103191_noctane_order_cryst_replicate.top
gmx mdrun 
# NVT-100K
gmx grompp -f nvt_100K.mdp -c topol.tpr -p ../103191_noctane_order_cryst_replicate.top -o new_topo.tpr
gmx mdrun -s new_topo.tpr -deffnm nvt_100K
# NOJUMP TRAJECTORY
echo 0 |gmx trjconv -f nvt_100K.xtc -s new_topo.tpr -pbc nojump -o nvt_100K_unwrap.xtc
