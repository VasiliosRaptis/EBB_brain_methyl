### make a .bz2 file with weights and metadata; to be used with fusion, see here: http://gusevlab.org/projects/fusion/

software=/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/software
EBB=/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation
name=EBB.BRAIN.METHYL.HERIT

mkdir -p $EBB/MWAS/$name
mkdir -p $EBB/MWAS/$name/$name
mkdir -p $EBB/MWAS/$name/scores/input
cd $EBB/MWAS/$name

## .list file
list_orig=$EBB/MWAS/$name.list
list_new=$name.list
cp $list_orig $list_new
sed -i "s#.*/##; s#^#$name/#" $list_new

## extract weights
while read -r wgtpath; do
  cp $wgtpath $EBB/MWAS/$name/$name
done < $list_orig

## .profile & .profile.err files
cp $EBB/MWAS/$name.profile* .

## .pos file - positions in hg37; P0: CpG pos -1; P1: CpG pos
Rscript $EBB/scripts/x_MWAS_04.5b_make_weigts_folder_helper.R 

## .variants file - all variant names
>$name.variants
for f in $EBB/MWAS/SCORES/input/* ; do
  sed '1d' $f >> $name.variants
done


## .sscore files
cd $EBB/MWAS/$name/scores/
cp -r $EBB/MWAS/SCORES/input/ .
# remove files with 0 variants
cd $EBB/MWAS/$name/scores/input
find . -name '*score.in' -exec wc -l {} + | awk '$1==1 {print $2}' > $EBB/MWAS/$name/temp.variants.to.remove
cd $EBB/MWAS/$name
while read -r badweight; do
  bname=$(basename $badweight)
  rm "scores/input/$bname" 
done < $EBB/MWAS/$name/temp.variants.to.remove
rm $EBB/MWAS/$name/temp.variants.to.remove

## create .tar.bz2 file
cd ../
tar -cvjSf $name.tar.bz2 $name
