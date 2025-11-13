#!bin/bash
nchains=50
natch=362
# First chain
idxs=1
idxe=353

echo "nmols: $nchains" >iPP_SC_40mon_50Ch_headtail.dat
echo "#ich idx-head idx-tail (indexes start at 0)" >>iPP_SC_40mon_50Ch_headtail.dat
for (( i=0; i<$nchains; i++ )) do
    echo "$i $idxs $idxe" >>iPP_SC_40mon_50Ch_headtail.dat
    idxs=`echo $idxs+$natch|bc -l`
    idxe=`echo $idxe+$natch|bc -l`
done
