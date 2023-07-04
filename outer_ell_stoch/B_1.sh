# Compute function $B^1$ related to outer map of elliptic problem, given by
#
#    B^1 = B_{in}^1 + B_{out}^1.

DATFILE=B_1.dat
RESFILE=B_1.res
ERRFILE=B_1.err

# Select only the first 117 lines, corresponding to energies H<-1.4854.
# The reason for this is because, for larger energies, g=0/PI is not a surface
# of section for the flow (either for po or hom.).

#lines=`wc -l ../outer_circ/omega_neg_SECg_br1.res | cut -d ' ' -f 1`
lines=116

echo $lines

cut -d ' ' -f 1,2 ../portbp/porbits.res >temp1   # H, T

cut -d ' ' -f 3 B_out_unst_br1.res >temp2 	# \Im(B^+)

# \Re(B_{in}), \Im(B_{in})
cut -d ' ' -f 2,3 ../inner_ell_stoch/inner_ell_stoch.res >temp3 

# \omega_-^j
cut -d ' ' -f 1 ../outer_circ/omega_neg_unst_br1.res >temp4 


echo "0.95387536e-3 0" > $DATFILE   # mu, UNSTABLE

# H, T, \Im(B^+), \Re(B_{in}), \Im(B_{in}), \omega_-^j
paste -d ' ' temp1 temp2 temp3 temp4 | sed -n '1,116p' >> $DATFILE
rm temp*

B_j <$DATFILE >$RESFILE

# now simply:
# - paste B_1 and B_2 onto B_12 using B_12.sh
# - plot the file B_12.res with gnuplot using B_12.plt
