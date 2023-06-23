FILEROOT=periodicorbits_del
DATFILE=$FILEROOT.dat
TMPFILE=$FILEROOT.tmp
RESFILE=$FILEROOT.res

cut -d ' ' -f 3-6 ../portbp/porbits.res >$DATFILE
./cardel <$DATFILE >$RESFILE
rm $DATFILE
