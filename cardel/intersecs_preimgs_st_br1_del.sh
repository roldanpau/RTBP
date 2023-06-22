FILEROOT=intersecs_preimgs_st_br1_del
DATFILE=$FILEROOT.dat
TMPFILE=$FILEROOT.tmp
RESFILE=$FILEROOT.res

cut -d ' ' -f 2-5 ../intersec/intersecs_st_br1.res >$DATFILE
./cardel <$DATFILE >$RESFILE
rm $DATFILE
