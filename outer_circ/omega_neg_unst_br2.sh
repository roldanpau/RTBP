# Compute integral $\omega_neg^f$ related to outer map for one of the
# homoclinic channels

NAMEROOT=omega_neg_unst_br2

datfile=$NAMEROOT.dat
resfile=$NAMEROOT.res
errfile=$NAMEROOT.err

echo "0.95387536e-3 SECg" >$datfile   # mu

cut -d ' ' -f 2 ../portbp/porbits.res > temp1    # period T

# z: homoclinic point z (in Delaunay)
cut -d ' ' -f 2-5 ../cardel/intersecs_unst_br2_del.res > temp2

# M: number of iterations to reach homoclinic point z from z_u
cut -d ' ' -f 2 ../approxint/approxints_unst_br2.res > temp3

# Select only the first 136 lines, corresponding to energies H<=-1.4494.
# The reason for this is because, for larger energies, the p.o. is not 
# in almost resonance anymore, and g=0/PI is not a surface of section.
paste -d ' ' temp1 temp2 temp3 | sed -n '1,136p' >>$datfile
rm temp1 temp2 temp3

./outer_circ_stoch < $datfile > $resfile 2> $errfile

# Warning: we don't want to remove datafile since outer_circ is executing in
# the background!
#rm $datfile
