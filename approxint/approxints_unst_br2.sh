cut -d ' ' -f 1-3 ../sec1sec2/sec1sec2.res > temp1.txt   	# H, p[2]
cut -d ' ' -f 4-5 ../hyper/hypers.res > temp2.txt      	# v_u[2]
cut -d ' ' -f 2 ../hyper/hypers.res > temp3.txt        	# rho_u

# Header: mu, k, unstable flag, branch flag, axis line $p_x=a$
echo "0.95387536e-3 2 0 0 0.0" > approxints.dat

paste -d ' ' temp1.txt temp2.txt temp3.txt >> approxints.dat
rm temp1.txt temp2.txt temp3.txt

./approxint <approxints.dat >approxints_unst_br2.res 2>approxints_unst_br2.err

rm approxints.dat
