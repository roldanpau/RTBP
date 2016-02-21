cut -d ' ' -f 1-3 ../sec1sec2/sec1sec2.res > temp1.txt   	# H, p[2]
cut -d ' ' -f 4-5 ../hyper/hypers.res > temp2.txt      	# v_u[2]
cut -d ' ' -f 2 ../hyper/hypers.res > temp3.txt        	# rho_u

# Header: mu, k, unstable flag, axis line $p_x=a$
echo "0.95387536e-3 2 0 0.0" > approxints.dat

paste -d ' ' temp1.txt temp2.txt temp3.txt | tac >> approxints.dat

rm temp1.txt temp2.txt temp3.txt

approxint <approxints.dat >approxints.res 2>approxints.err
tac approxints.res >approxints.res.tmp
mv approxints.res.tmp approxints.res

rm approxints.dat
