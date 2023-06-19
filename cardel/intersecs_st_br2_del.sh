FILEROOT=intersecs_st_br2_del
DATFILE=$FILEROOT.dat
TMPFILE=$FILEROOT.tmp
RESFILE=$FILEROOT.res

echo "0.95387536e-3" >$DATFILE
cut -d ' ' -f 1,5 ../intersec/intersecs_st_br2.res >$TMPFILE
awk '{print $1, $2, 0}' $TMPFILE >>$DATFILE
./cardels_2d <$DATFILE >$RESFILE
rm $TMPFILE $DATFILE
