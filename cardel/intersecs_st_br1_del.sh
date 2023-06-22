FILEROOT=intersecs_st_br1_del
DATFILE=$FILEROOT.dat
TMPFILE=$FILEROOT.tmp
RESFILE=$FILEROOT.res

cut -d ' ' -f 7-10 ../intersec/intersecs_st_br1.res >$DATFILE
./cardel <$DATFILE >$RESFILE
rm $DATFILE
