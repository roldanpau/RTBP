NAMEROOT=intersecs_preimgs_unst_br1_SECg2

datfile=$NAMEROOT.dat
resfile=$NAMEROOT.res
errfile=$NAMEROOT.err

echo "0.95387536e-3 SECg2 1" > $datfile
cat ../cardel/intersecs_preimgs_unst_br1_del.res >> $datfile
./prtbpdel < $datfile > $resfile
rm $datfile
