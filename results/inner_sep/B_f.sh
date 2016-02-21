# Compute function $B^{f,+}$ related to outer map of elliptic problem, given
# by
#
#    B^f = B_{out}^f + B_{in}^f e^{i\mu\omega_{out}^f},
#
# where B_{out}^f = i 2\mu\Im(B^+), 	
# and \omega_{out}^f = 2\omega_+^f.

# remember that, since Marcel changed signs in the last version of the paper,
# our computed value of \omega_+ is now the negative of Marcel's \omega_+.


lines=`wc -l omega_pos_f.res | cut -d ' ' -f 1`

echo $lines

# Select only the first $lines lines, corresponding to energy levels H=-1.56
# -> -1.827
cut -d ' ' -f 1 ../../portbp/porbits.res | tac | head -n $lines >temp1 	# H

cut -d ' ' -f 2 B_out_f.res >temp2 	# \Im(B^+)

# Select only the first $lines lines, corresponding to energy levels H=-1.56
# -> -1.827
# \Re(B_{in}^f), \Im(B_{in}^f)
cut -d ' ' -f 4,5 ../inner_ell.res | head -n $lines >temp3 

# H, \Im(B^+), \Re(B_{in}^f), \Im(B_{in}^f), \omega_pos
paste -d ' ' temp1 temp2 temp3 omega_pos_f.res > B_f.dat 	
rm temp*

B_f <B_f.dat >B_f.res

# now simply plot the file B_f.res with gnuplot using B_f.plt
