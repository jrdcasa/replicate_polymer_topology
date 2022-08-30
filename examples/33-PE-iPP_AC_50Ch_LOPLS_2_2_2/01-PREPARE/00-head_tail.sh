#!bin/bash
nchains=50
nchA=25
nchB=25
natchA=242
# First kind chain
idxsA=1
idxeA=236

echo "nmols: $nchains" >iPP-PE_40mon_25-25Ch_headtail.dat
echo "#ich idx-head idx-tail (indexes start at 0)" >>iPP-PE_40mon_25-25Ch_headtail.dat
for (( i=0; i<$nchA; i++ )) do
    echo "$i $idxsA $idxeA" >>iPP-PE_40mon_25-25Ch_headtail.dat
    idxsA=`echo $idxsA+$natchA|bc -l`
    idxeA=`echo $idxeA+$natchA|bc -l`
done
# Second kind chain
natchB=362
idxsB=6051
idxeB=6403
for (( i=$nchA; i<$nchains; i++ )) do
    echo "$i $idxsB $idxeB" >>iPP-PE_40mon_25-25Ch_headtail.dat
    idxsB=`echo $idxsB+$natchB|bc -l`
    idxeB=`echo $idxeB+$natchB|bc -l`
done
