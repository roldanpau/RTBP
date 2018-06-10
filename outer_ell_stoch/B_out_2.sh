# Compute function $B_{out}^2$ associated to outer map of elliptic problem.
#
# Actually, we only compute the integral called $B^+$ in my notes. This is
# enough, because 
#    B_{out}^j = -\mu B^+ + \mu C^+$,
# where
#    re(C^+)=re(B^+), and im(C^+)=-im(B^+).

DATFILE=B_out_2.dat
RESFILE=B_out_2.res
ERRFILE=B_out_2.err

lines=`wc -l ../outer_circ/omega_neg_SECg2_br1.res | cut -d ' ' -f 1`

echo $lines

# Select only the first $lines lines, corresponding to energy levels H= -> 
cut -d ' ' -f 1 ../portbp/porbits.res | head -n $lines >temp1 	# H

# Select only the first $lines lines, corresponding to energy levels H= ->
cut -d ' ' -f 1-4 ../prtbp_del_car/prtbp_del_cars_SECg2.res | head -n $lines \
    > temp2  # p

# z_u: preimage of homoclinic point z
cut -d ' ' -f 13-16 ../intersec_del_car/intersecs_unst_SECg2_br1.res \
   | head -n $lines > temp3

# z_s: preimage of homoclinic point z
# We could use 
#   cut -d ' ' -f 13-16 ../intersec_del_car/intersecs_st_SECg_br1.res > temp3
# but we decide to simply use the symmetry z_s = z_u, but reversing the sign of
# l and g. (This is done in the program).

# M: number of iterations to reach homoclinic point z from z_u
cut -d ' ' -f 2 ../approxint_del_car/approxints_unst_SECg2_br1.res \
    | head -n $lines > temp4

echo "0.95387536e-3" > $DATFILE

# H, p, z_u, \omega_pos, N
paste -d ' ' temp1 temp2 temp3 ../outer_circ/omega_neg_SECg2_br1.res temp4 \
    >> $DATFILE
rm temp*

nohup ./outer_ell_stoch <$DATFILE >$RESFILE 2>$ERRFILE &
