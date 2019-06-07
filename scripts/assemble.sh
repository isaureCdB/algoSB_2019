nfrag=$1
cutoff=$2
npose=$3
meanrank=$4
maxchains=$5

if [ ! -s pairs_${cutoff}A_${npose}poses.json ];then

  echo ""
  echo '**********************************************************'
  echo 'Extract docking coordinates'
  echo '**********************************************************'
  list_cg=partner2-ensemble.list
  list_aa=partner2-ensemble-aa.list

  template_cg=`head -n 1 $list_cg`
  template_aa=`head -n 1 $list_aa`

  #convert results into binary coordinates
  python2 $ATTRACTDIR/dump_coordinates.py result.dat $template_cg AAA-preatoms.npy `cat $ALGOSB/exo3/AAA.preatoms` --ens 2 $list_cg
  python2 $ATTRACTDIR/dump_coordinates.py result.dat $template_cg AAA-postatoms.npy `cat $ALGOSB/exo3/AAA.postatoms` --ens 2 $list_cg
  python2 $ATTRACTDIR/dump_coordinates.py result.dat $template_aa AAA.npy -1 --ens 2 $list_aa

  for frag in `seq $nfrag`; do
    ln -s AAA-preatoms.npy frag$frag-preatoms.npy
    ln -s AAA-postatoms.npy frag$frag-postatoms.npy
  done

  echo ""
  echo '**********************************************************'
  echo 'Calculate 2-fragments connections'
  echo '**********************************************************'
  #name="${cutoff}A_${npose}poses"
  echo ${cutoff}A_${npose}poses
  $ALGOSB/scripts/connect.py 2 $cutoff $npose 200 AAA-preatoms.npy AAA-preatoms.npy AAA-postatoms.npy AAA-postatoms.npy  > pairs_${cutoff}A_${npose}poses.json
fi

a=`cat check_${cutoff}A_${npose}poses`
if [ $a -lt 1 ];then exit ; fi

echo ""
echo '**********************************************************'
echo 'Propagate into a $nfrags-fragments connection graph'
echo '**********************************************************'
if [ ! -s graph_${nfrag}frag_${cutoff}A_${npose}poses.json ];then
  $ALGOSB/scripts/connect-homo.py pairs_${cutoff}A_${npose}poses.json $nfrag > graph_${nfrag}frag_${cutoff}A_${npose}poses.json
fi

a="frag[1-"$nfrag"]-preatoms.npy"
b="frag[1-"$nfrag"]-postatoms.npy"

name=chains_${nfrag}frag_${cutoff}A_${npose}poses

if [ ! -s $name.list ];then
  echo ""
  echo '**********************************************************'
  echo 'Sample chains'
  echo '**********************************************************'
  $ALGOSB/scripts/make_chains.py graph_${nfrag}frag_${cutoff}A_${npose}poses.json \
      --meanrank $meanrank --preatoms $a --postatoms $b  --maxchains $maxchains \
      |awk 'NR>1{for (i=4;i<=NF;i++) printf("%i ",$i); printf("\n")}' \
      > $name.list
fi

if [ ! -s $name.list ];then
  echo "try increasing npose or meanrank or cutoff, or decrease nfrag"
  exit
fi

echo ""
echo '**********************************************************'
echo 'Build chains structure'
echo '**********************************************************'

name=chains_${nfrag}frag_${cutoff}A_${npose}poses

# merge the overlapping nucleotides:
  # nat = nb of atom per nucleotide
  # motifs = index of motif(sequence) for each fragment
  # npy = coordinates of each motif

if [ ! -s $name.npy ];then
  $ALGOSB/scripts/chain2rna.py $name.list \
    --nfrag $nfrag \
    --nat `cat $ALGOSB/exo3/nat.list` \
    --motifs `cat $ALGOSB/exo3/motifs.list` \
    --npy AAA.npy \
    --outp $name.npy
fi

#compute the rmsd of the chains
$ALGOSB/scripts/rmsdnpy.py $name.npy $ALGOSB/exo3/rna_b.pdb > $name.rmsd

echo "RMSD rank"
awk '{print $2, NR}' $name.rmsd |sort -nk1|head -n 5

#select chains with RMSD < 5
awk '$2<5{print NR}' $name.rmsd > inf5A
a=`cat inf5A|wc -l`
if [ $a -lt 1 ];then exit ; fi

echo ""
echo "**********************************************************"
echo " prepare chains with RMSD < 5 for visualisation"
echo "**********************************************************"
$ALGOSB/scripts/select-npy.py $name.npy $name-inf5A.npy --struct `cat inf5A`

#convert binary coordinates into readable coordinates
nnucl=$(($nfrag+2))
template_rna=$ALGOSB/exo3/template/A${nnucl}.pdb
$ALGOSB/scripts/npy2pdb.py $name-inf5A.npy $template_rna | \
awk 'substr($3,1,1)!="H"' > $name-inf5A.pdb
