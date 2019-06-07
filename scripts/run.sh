name=$1

tar xzf  $name.tgz
cd $name/
sed -i 's/receptor\-heavy\.pdb\ ligand\-heavy\.pdb/\/dev\/null\/\ ligand\-heavy\.pdb/' $name.sh
./$name.sh

if [ -s results.irmsd ]; then
  echo "************************************"
  echo "  top iRMSD values"
  echo "************************************"
  echo "iRMSD rank"
  awk '{print $2, NR}' result.irmsd|sort -nk1|head

  $ALGOSB/scripts/plot-rmsd-rank.sh result.irmsd

fi
rm receptorgrid.gridheader
rm out_$name.dat out_$name.score out_$name-scored.dat out_$name-sorted.dat
rm systsearch*dat
