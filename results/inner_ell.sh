# Compute integral A^+(I) related to inner map of elliptic problem.
# Moreover, compute integral B_ín^{f,+}(I) related to outer map of elliptic
# problem.

# We output lines with INCREASING energy H=-1.56 -> -2.04

echo "0.95387536e-3" >inner_ell.dat	# mu

# reverse lines, to get INCREASING energy values
cut -d ' ' -f 1 ../portbp/porbits.res | tac > temp1 	# H
cut -d ' ' -f 1-4 prtbpdel_p3.res | tac > temp2 	# periodic point p3
cut -d ' ' -f 1-4 prtbpdel_p4.res | tac > temp3 	# periodic point p4
paste -d ' ' temp1 temp2 temp3 >> inner_ell.dat
rm temp1 temp2 temp3 

nohup inner_ell < inner_ell.dat > inner_ell.res 2> inner_ell.err &
