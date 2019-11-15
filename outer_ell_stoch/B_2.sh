# Compute function $B^2$ related to outer map of elliptic problem, given by
#
#    B^2 = B_{in}^2 + B_{out}^2.

DATFILE=B_2.dat
RESFILE=B_2.res
ERRFILE=B_2.err

lines=`wc -l ../outer_circ/omega_neg_SECg2_br1.res | cut -d ' ' -f 1`

echo $lines

# Select only the first $lines lines, corresponding to energy levels H= -> 
cut -d ' ' -f 1,2 ../portbp/porbits.res | head -n $lines >temp1   # H, T

cut -d ' ' -f 3 B_out_2.res >temp2 	# \Im(B^+)

# Select only the first $lines lines, corresponding to energy levels H= ->

# \Re(B_{in}), \Im(B_{in})
cut -d ' ' -f 2,3 ../inner_ell_stoch/inner_ell_stoch.res | head -n $lines >temp3 

# \omega_-^2
cut -d ' ' -f 1 ../outer_circ/omega_neg_SECg2_br1.res >temp4 


# H, T, \Im(B^+), \Re(B_{in}), \Im(B_{in}), \omega_-^j
paste -d ' ' temp1 temp2 temp3 temp4 > $DATFILE
rm temp*

B_j <$DATFILE >$RESFILE

# now simply:
# - paste B_1 and B_2 onto B_12 using B_12.sh
# - plot the file B_12.res with gnuplot using B_12.plt
