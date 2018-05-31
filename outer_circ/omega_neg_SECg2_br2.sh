# Compute integral $\omega_neg^j$ related to outer map for the POINCARE SECTION
# {g=\pi}.

datfile=omega_neg_SECg2_br2.dat
resfile=omega_neg_SECg2_br2.res
errfile=omega_neg_SECg2_br2.err

echo "0.95387536e-3 SECg2" >$datfile   # mu

cut -d ' ' -f 2 ../portbp/porbits.res > temp1    # period T

# z_u, z_u_car: preimage of homoclinic point z
cut -d ' ' -f 13-20 ../intersec_del_car/intersecs_unst_SECg2_br2.res > temp2

# M: number of iterations to reach homoclinic point z from z_u
cut -d ' ' -f 2 ../approxint_del_car/approxints_unst_SECg2_br2.res > temp3

paste -d ' ' temp1 temp2 temp3 >>$datfile
rm temp1 temp2 temp3

./outer_circ_stoch < $datfile > $resfile 2> $errfile
# Warning: we don't want to remove datafile since outer_circ is executing in
# the background!
#rm $datfile
