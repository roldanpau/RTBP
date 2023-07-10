# Compute first order of the variance, $\sigma_0^2$, given by

fileroot=variance
DATFILE=$fileroot.dat
RESFILE=$fileroot.res
ERRFILE=$fileroot.err

# Select only the lines corresponding to energies -1.551 <= H <= -1.4894.
# This corresponds to interval I_1 in the paper (Table 2).

cut -d ' ' -f 1 B_2_I1.dat >temp1   # H
cut -d ' ' -f 5,6 B_2_I1.dat >temp2   # Re(B1), Im(B1)
cut -d ' ' -f 5,6 B_3_I1.dat >temp3   # Re(B2), Im(B2)

# H, alpha_neg_1, alpha_pos_2, Re(B1), Im(B1), Re(B2), Im(B2)
paste -d ' ' temp1 alpha_neg_2_I1.dat alpha_pos_3_I1.dat temp2 temp3 > $DATFILE
rm temp*

./variance <$DATFILE >$RESFILE
