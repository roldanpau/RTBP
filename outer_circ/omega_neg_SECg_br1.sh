# Compute integral $\omega_neg^f$ related to outer map for the POINCARE SECTION
# {g=0}.

datfile=omega_neg_SECg_br1.dat
resfile=omega_neg_SECg_br1.res
errfile=omega_neg_SECg_br1.err

echo "0.95387536e-3 SECg" >$datfile   # mu

cut -d ' ' -f 2 ../portbp/porbits.res > temp1    # period T

# z_u, z_u_car: preimage of homoclinic point z
cut -d ' ' -f 13-20 ../intersec_del_car/intersecs_unst_SECg_br1.res > temp2

# M: number of iterations to reach homoclinic point z from z_u
cut -d ' ' -f 2 ../approxint_del_car/approxints_unst_SECg_br1.res > temp3

# Choose one: 
# 1. Select only the first 116 lines, corresponding to energies H<=-1.4894
# 2. Select only the first 136 lines, corresponding to energies H<=-1.4494.
# The reason for this is because, for larger energies, the p.o. is not 
# in almost resonance anymore, and g=0/PI is not a surface of section.
paste -d ' ' temp1 temp2 temp3 | sed -n '1,136p' >>$datfile
rm temp1 temp2 temp3

./outer_circ_stoch < $datfile > $resfile 2> $errfile

# Warning: we don't want to remove datafile since outer_circ is executing in
# the background!
#rm $datfile
