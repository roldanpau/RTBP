FILEROOT=splittings_st_br1
DATFILE=$FILEROOT.dat
RESFILE=$FILEROOT.res
ERRFILE=$FILEROOT.err

cut -d ' ' -f 1 ../sec1sec2/sec1sec2_2d.res  > temp1.txt  # H,
cut -d ' ' -f 6-7 ../hyper/hypers.res  > temp2.txt      	# v_s[2]
cut -d ' ' -f 2 ../approxint/approxints_st_br1.res  > temp3.txt   # n
cut -d ' ' -f 2-3 ../intersec/intersecs_st_br1.res > temp4.txt   # p_s[2]

# Header: mu, stable flag, RIGHT branch
echo "0.95387536e-3 1 1" > $DATFILE

paste -d ' ' temp1.txt temp2.txt temp3.txt temp4.txt >> $DATFILE

rm temp1.txt temp2.txt temp3.txt temp4.txt

./splitting <$DATFILE >$RESFILE 2>$ERRFILE

#rm $DATFILE
