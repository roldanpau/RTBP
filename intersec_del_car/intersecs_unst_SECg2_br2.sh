datfile=intersecs_unst_SECg2_br2.dat
resfile=intersecs_unst_SECg2_br2.res
errfile=intersecs_unst_SECg2_br2.err
tmpfile=intersecs_unst_SECg2_br2.res.tmp

cut -d ' ' -f 1-3 ../sec1sec2/sec1sec2_2d.res > temp1.txt   	# H, p[2]
cut -d ' ' -f 4-5 ../hyper/hypers.res > temp2.txt      	# v_u[2]
cut -d ' ' -f 2 ../hyper/hypers.res > temp3.txt        	# rho_u

cut -d ' ' -f 2-4 ../approxint_del_car/approxints_unst_SECg2_br2.res \
    > temp4.txt    # n, h1, h2


# Header: mu, sec, branch, unstable flag, axis line $g=a$
echo "0.95387536e-3 SECg2 LEFT 0 0.0" > $datfile

paste -d ' ' temp1.txt temp2.txt temp3.txt temp4.txt | tac >> $datfile

rm temp1.txt temp2.txt temp3.txt temp4.txt

./intersec_del_car <$datfile >$resfile 2>$errfile
tac $resfile >$tmpfile
mv $tmpfile $resfile

#rm $datfile
