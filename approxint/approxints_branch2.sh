cut -d ' ' -f 1-3 ../sec1sec2/sec1sec2.res > temp1.txt   	# H, p[2]
cut -d ' ' -f 6-7 ../hyper/hypers.res > temp2.txt      	# v_s[2]
cut -d ' ' -f 3 ../hyper/hypers.res > temp3.txt        	# rho_s

# Header: mu, k, stable flag, axis line $p_x=a$
echo "0.95387536e-3 2 1 0.0" > approxints_branch2.dat

paste -d ' ' temp1.txt temp2.txt temp3.txt | tac >> approxints_branch2.dat

rm temp1.txt temp2.txt temp3.txt

approxint <approxints_branch2.dat >approxints_branch2.res 2>approxints_branch2.err
tac approxints_branch2.res >approxints_branch2.res.tmp
mv approxints_branch2.res.tmp approxints_branch2.res

rm approxints_branch2.dat
