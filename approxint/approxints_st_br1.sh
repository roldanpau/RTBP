NAMEROOT=approxints_st_br1
DATFILE=$NAMEROOT.dat
RESFILE=$NAMEROOT.res
ERRFILE=$NAMEROOT.err

cut -d ' ' -f 1-3 ../sec1sec2/sec1sec2_2d.res > temp1.txt   	# H, p[2]
cut -d ' ' -f 6-7 ../hyper/hypers.res > temp2.txt      	# v_s[2]
cut -d ' ' -f 3 ../hyper/hypers.res > temp3.txt        	# rho_s

# Header: mu, k, stable flag, branch flag, axis line $p_x=a$
echo "0.95387536e-3 2 1 1 0.0" > $DATFILE

paste -d ' ' temp1.txt temp2.txt temp3.txt >> $DATFILE
rm temp1.txt temp2.txt temp3.txt

./approxint <$DATFILE >$RESFILE 2>$ERRFILE

#rm $DATFILE
