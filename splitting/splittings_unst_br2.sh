FILEROOT=splittings_unst_br2
DATFILE=$FILEROOT.dat
RESFILE=$FILEROOT.res
ERRFILE=$FILEROOT.err

cut -d ' ' -f 1 ../sec1sec2/sec1sec2_2d.res  > temp1.txt  # H,
cut -d ' ' -f 4-5 ../hyper/hypers.res  > temp2.txt      	# v_u[2]
cut -d ' ' -f 2 ../approxint/approxints_unst_br2.res  > temp3.txt   # n
cut -d ' ' -f 2-3 ../intersec/intersecs_unst_br2.res > temp4.txt   # p_u[2]

# Header: mu, unstable flag, LEFT branch
echo "0.95387536e-3 0 0" > $DATFILE

paste -d ' ' temp1.txt temp2.txt temp3.txt temp4.txt >> $DATFILE

rm temp1.txt temp2.txt temp3.txt temp4.txt

./splitting <$DATFILE >$RESFILE 2>$ERRFILE

#rm $DATFILE
