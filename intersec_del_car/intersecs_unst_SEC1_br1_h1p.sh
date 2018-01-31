datfile=intersecs_unst_SEC1_br1_h1p.dat
resfile=intersecs_unst_SEC1_br1_h1p.res
errfile=intersecs_unst_SEC1_br1_h1p.err
tmpfile=intersecs_unst_SEC1_br1_h1p.res.tmp

cut -d ' ' -f 1-3 ../sec1sec2/sec1sec2.res > temp1.txt   	# H, p[2]
cut -d ' ' -f 4-5 ../hyper/hypers.res > temp2.txt      	# v_u[2]
cut -d ' ' -f 2 ../hyper/hypers.res > temp3.txt        	# rho_u

cut -d ' ' -f 2-4 ../approxint_del_car/approxints_unst_SEC1_br1.res > temp4.txt    # n, h1, h2


# Header: mu, unstable flag, axis line $g=a$
echo "0.95387536e-3 SEC1 0 3.14160265358979323844" > $datfile

paste -d ' ' temp1.txt temp2.txt temp3.txt temp4.txt |tac >> $datfile

rm temp1.txt temp2.txt temp3.txt temp4.txt

./intersec_del_car <$datfile >$resfile 2>$errfile
tac $resfile >$tmpfile
mv $tmpfile $resfile

rm $datfile
