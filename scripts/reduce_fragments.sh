for i in `seq 8`; do
  python $ATTRACTTOOLS/reduce.py frag$i.pdb --rna
done
