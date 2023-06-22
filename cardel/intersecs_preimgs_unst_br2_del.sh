FILEROOT=intersecs_preimgs_unst_br2_del
DATFILE=$FILEROOT.dat
TMPFILE=$FILEROOT.tmp
RESFILE=$FILEROOT.res

cut -d ' ' -f 2-5 ../intersec/intersecs_unst_br2.res >$DATFILE
./cardel <$DATFILE >$RESFILE
rm $DATFILE
