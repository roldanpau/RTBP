cut -d ' ' -f 1-3 ../sec1sec2/sec1sec2_notall.res > temp1.txt   	# H, p[2]
cut -d ' ' -f 6-7 ../hyper/hypers_notall.res > temp2.txt      	# v_s[2]
cut -d ' ' -f 3 ../hyper/hypers_notall.res > temp3.txt        	# rho_s

cut -d ' ' -f 2-4 ../approxint/approxints_branch2.res > temp4.txt    # n, h1, h2


# Header: mu, stable flag, axis line $p_x=a$
echo "0.95387536e-3 1 0.0" > intersecs_branch2.dat

paste -d ' ' temp1.txt temp2.txt temp3.txt temp4.txt |tac >> intersecs_branch2.dat

rm temp1.txt temp2.txt temp3.txt temp4.txt

intersec <intersecs_branch2.dat >intersecs_branch2.res
tac intersecs_branch2.res >intersecs_branch2.res.tmp
mv intersecs_branch2.res.tmp intersecs_branch2.res

rm intersecs_branch2.dat
