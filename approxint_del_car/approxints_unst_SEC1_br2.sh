datfile=approxints_unst_SEC1_br2.dat
resfile=approxints_unst_SEC1_br2.res
errfile=approxints_unst_SEC1_br2.err

cut -d ' ' -f 1-3 ../sec1sec2/sec1sec2.res > temp1.txt   	# H, p[2]
cut -d ' ' -f 4-5 ../hyper/hypers.res > temp2.txt      	# v_u[2]
cut -d ' ' -f 2 ../hyper/hypers.res > temp3.txt        	# rho_u

# Header: mu, sec, k, unstable flag, branch flag, axis line $p_x=a$
echo "0.95387536e-3 SEC1 2 0 0 3.14159265358979323844" > $datfile

paste -d ' ' temp1.txt temp2.txt temp3.txt >> $datfile
rm temp1.txt temp2.txt temp3.txt

./approxint_del_car <$datfile >$resfile 2>$errfile

rm $datfile
