cut -d ' ' -f 1-3 ../sec1sec2/sec1sec2.res > temp1.txt   	# H, p[2]
cut -d ' ' -f 4-5 ../hyper/hypers.res > temp2.txt      	# v_u[2]
cut -d ' ' -f 2 ../hyper/hypers.res > temp3.txt        	# rho_u

cut -d ' ' -f 2-4 ../approxint_del_car/approxints_unst_br1.res > temp4.txt    # n, h1, h2


# Header: mu, unstable flag, axis line $g=a$
echo "0.95387536e-3 0 0.0" > intersecs_unst_br1.dat

paste -d ' ' temp1.txt temp2.txt temp3.txt temp4.txt |tac >> intersecs_unst_br1.dat

rm temp1.txt temp2.txt temp3.txt temp4.txt

./intersec_del_car <intersecs_unst_br1.dat >intersecs_unst_br1.res 2>intersecs_unst_br1.err
tac intersecs_unst_br1.res >intersecs_unst_br1.res.tmp
mv intersecs_unst_br1.res.tmp intersecs_unst_br1.res

rm intersecs_unst_br1.dat
