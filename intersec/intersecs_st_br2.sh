FILEROOT=intersecs_st_br2
DATFILE=$FILEROOT.dat
RESFILE=$FILEROOT.res
ERRFILE=$FILEROOT.err

cut -d ' ' -f 1-3 ../sec1sec2/sec1sec2_2d.res > temp1.txt   	# H, p[2]
cut -d ' ' -f 6-7 ../hyper/hypers.res > temp2.txt      	# v_s[2]
cut -d ' ' -f 3 ../hyper/hypers.res > temp3.txt        	# rho_s

cut -d ' ' -f 2-4 ../approxint/approxints_st_br2.res > temp4.txt    # n, h1, h2


# Header: mu, stable flag, axis line $p_x=a$
echo "0.95387536e-3 1 0.0" > $DATFILE

paste -d ' ' temp1.txt temp2.txt temp3.txt temp4.txt >> $DATFILE

rm temp1.txt temp2.txt temp3.txt temp4.txt

./intersec <$DATFILE >$RESFILE 2>$ERRFILE

#rm $DATFILE
