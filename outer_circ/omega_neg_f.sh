# Compute integral $\omega_neg^f$ related to outer map.

echo "0.95387536e-3" >omega_neg_f.dat   # mu

cut -d ' ' -f 2 ../portbp/porbits_notall.res > temp1    # period T

# z_u: preimage of homoclinic point z
cut -d ' ' -f 1-4 ../prtbp_del/prtbpdel_zu.res > temp2

# M: number of iterations to reach homoclinic point z from z_u
cut -d ' ' -f 2 ../approxint/approxints.res > temp3

paste -d ' ' temp1 temp2 temp3 >>omega_neg_f.dat
rm temp1 temp2 temp3

./outer_circ < omega_neg_f.dat > omega_neg_f.res 2> omega_neg_f.err &
# Warning: we don't want to remove datafile since outer_circ is executing in
# the background!
#rm omega_neg_f.dat
