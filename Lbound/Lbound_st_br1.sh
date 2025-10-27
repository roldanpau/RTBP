# Compute \f$L\f$ bound for one of the homoclinic channels

NAMEROOT=Lbound_st_br1

datfile=$NAMEROOT.dat
resfile=$NAMEROOT.res
errfile=$NAMEROOT.err

echo "0.95387536e-3 1" >$datfile   # mu, STABLE

cut -d ' ' -f 2 ../portbp/porbits.res > temp1    # period T

# zu: preimage of homoclinic point z (in Cartesian)
cut -d ' ' -f 2-5 ../intersec/intersecs_st_br1.res > temp2

# t: integration time to reach homoclinic point z from z_u
cut -d ' ' -f 6 ../intersec/intersecs_st_br1.res > temp3

# Select only the first 117 lines, corresponding to energies H<-1.4854.
# The reason for this is because, for larger energies, g=0/PI is not a surface
# of section for the flow (either for po or hom.).
paste -d ' ' temp1 temp2 temp3 | sed -n '1,117p' >>$datfile
rm temp1 temp2 temp3

./Lbound < $datfile > $resfile.tmp 2> $errfile

# Add column with energy H to results file
cut -d ' ' -f 1 ../portbp/porbits.res > temp1    # energy H
paste -d ' ' temp1 $resfile.tmp | sed -n '1,117p' > $resfile
rm temp1
rm $resfile.tmp

rm $datfile
