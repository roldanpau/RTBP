# Compute integral A_ín(I) related to outer map of elliptic problem.

echo "0.95387536e-3" >inner_ell_stoch.dat	# mu

cut -d ' ' -f 1 ../portbp/porbits.res > temp1 	# H
cut -d ' ' -f 1-4 ../prtbp_del_car/prtbp_del_cars_SECg.res > temp2 	# periodic point p
paste -d ' ' temp1 temp2 >> inner_ell_stoch.dat
rm temp1 temp2 

nohup inner_ell_stoch < inner_ell_stoch.dat > inner_ell_stoch.res \
    2> inner_ell_stoch.err &
