cut -d ' ' -f 1-3 ../sec1sec2/sec1sec2_notall.res > temp1.txt   	# H, p[2]
cut -d ' ' -f 4-5 ../hyper/hypers_notall.res > temp2.txt      	# v_u[2]
cut -d ' ' -f 2 ../hyper/hypers_notall.res > temp3.txt        	# rho_u

cut -d ' ' -f 2-4 ../approxint/approxints.res > temp4.txt    # n, h1, h2


# Header: mu, unstable flag, axis line $p_x=a$
echo "0.95387536e-3 0 0.0" > intersecs.dat

paste -d ' ' temp1.txt temp2.txt temp3.txt temp4.txt |tac >> intersecs.dat

rm temp1.txt temp2.txt temp3.txt temp4.txt

intersec <intersecs.dat >intersecs.res
tac intersecs.res >intersecs.res.tmp
mv intersecs.res.tmp intersecs.res

rm intersecs.dat
