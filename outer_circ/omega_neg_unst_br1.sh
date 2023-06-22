# Compute integral $\omega_neg^f$ related to outer map for one of the
# homoclinic channels

# REMEMBER TO CHANGE omega_pos_stoch -> omega_neg_stoch IN FILE
# outer_circ_stoch.c BEFORE CALLING THIS FILE!!!

NAMEROOT=omega_neg_unst_br1

datfile=$NAMEROOT.dat
resfile=$NAMEROOT.res
errfile=$NAMEROOT.err

echo "0.95387536e-3" >$datfile   # mu

cut -d ' ' -f 2 ../portbp/porbits.res > temp1    # period T

# zu: preimage of homoclinic point z (in Cartesian)
cut -d ' ' -f 2-5 ../intersec/intersecs_unst_br1.res > temp2

# t: integration time to reach homoclinic point z from z_u
cut -d ' ' -f 6 ../intersec/intersecs_unst_br1.res > temp3

# Select only the first 136 lines, corresponding to energies H<=-1.4494.
# The reason for this is because, for larger energies, the p.o. is not 
# in almost resonance anymore, and g=0/PI is not a surface of section.
paste -d ' ' temp1 temp2 temp3 | sed -n '1,136p' >>$datfile
rm temp1 temp2 temp3

./outer_circ_stoch < $datfile > $resfile 2> $errfile

# Warning: we don't want to remove datafile since outer_circ is executing in
# the background!
#rm $datfile
