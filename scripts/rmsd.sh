
bound=$ALGOSB/exo3/rna_b-aa.pdb

egrep "RA B   [123]" $bound > frag1_b.pdb
egrep "RA B   [234]" $bound > frag2_b.pdb
egrep "RA B   [345]" $bound > frag3_b.pdb
egrep "RA B   [456]" $bound > frag4_b.pdb
egrep "RA B   [567]" $bound > frag5_b.pdb
egrep "RA B   [678]" $bound > frag6_b.pdb

for i in `seq 6`; do

    python2 $ATTRACTDIR/lrmsd.py result.dat frag$i\_b.pdb frag$i\_b.pdb \
        --ens 2 partner2-ensemble-aa.list --allatoms  > frag$i.rmsd

    echo "**********************************"
    echo "  frag $i  "
    echo "**********************************"
    echo "RMSD rank"
    awk '{print $2, NR}' frag$i.rmsd |sort -nk1 > frag$i.rmsd-sorted
    head -n 5 frag$i.rmsd-sorted

done
