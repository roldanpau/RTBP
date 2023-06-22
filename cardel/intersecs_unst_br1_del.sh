FILEROOT=intersecs_unst_br1_del
DATFILE=$FILEROOT.dat
TMPFILE=$FILEROOT.tmp
RESFILE=$FILEROOT.res

cut -d ' ' -f 7-10 ../intersec/intersecs_unst_br1.res >$DATFILE
./cardel <$DATFILE >$RESFILE
rm $DATFILE
