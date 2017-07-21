# Compute integral $\omega_neg^f$ related to outer map for the FIRST POINCARE
# SECTION {l=0}.

datfile=omega_neg_f_SEC1.dat
resfile=omega_neg_f_SEC1.res
errfile=omega_neg_f_SEC1.err

echo "0.95387536e-3" >$datfile   # mu

cut -d ' ' -f 2 ../portbp/porbits.res > temp1    # period T

# z_u: preimage of homoclinic point z
cut -d ' ' -f 13-16 ../intersec_del_car/intersecs_unst_SEC1_br1.res > temp2

# M: number of iterations to reach homoclinic point z from z_u
cut -d ' ' -f 2 ../approxint_del_car/approxints_unst_SEC1_br1.res > temp3

paste -d ' ' temp1 temp2 temp3 >>$datfile
rm temp1 temp2 temp3

./outer_circ < $datfile > $resfile 2> $errfile
# Warning: we don't want to remove datafile since outer_circ is executing in
# the background!
#rm $datfile
