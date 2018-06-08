datfile=approxints_st_SECg2_br1.dat
resfile=approxints_st_SECg2_br1.res
errfile=approxints_st_SECg2_br1.err

cut -d ' ' -f 1-3 ../sec1sec2/sec1sec2_2d.res > temp1.txt   	# H, p[2]
cut -d ' ' -f 6-7 ../hyper/hypers.res > temp2.txt      	# v_s[2]
cut -d ' ' -f 3 ../hyper/hypers.res > temp3.txt        	# rho_s

# Header: mu, sec, k, stable flag, branch flag, axis line $g=a$
# Note: k refers to iterates on the Euclidean section, NOT on SECg.
# Note: the axis line $g=a$ is actually NOT used for section SECg.
echo "0.95387536e-3 SECg2 2 1 1 0.0" > $datfile

paste -d ' ' temp1.txt temp2.txt temp3.txt >> $datfile
rm temp1.txt temp2.txt temp3.txt

./approxint_del_car <$datfile >$resfile 2>$errfile

#rm $datfile
