# Compute function $B^{b,+}$ related to outer map of elliptic problem, given
# by
#
#    B^b = B_{in}^b + B_{out}^b e^{i\mu\omega_{in}^b},
#
# where B_{out}^b = i 2\mu\Im(B^+).

# remember that we have changed signs in the code so that our computed value
# of \omega_{in}^b is the same as Marcel's.


lines=`wc -l omega_b.dat | cut -d ' ' -f 1`

echo $lines

# Select only the first $lines lines, corresponding to energy levels H=-1.56
# -> -1.827
cut -d ' ' -f 1 ../portbp/porbits.res | tac | head -n $lines >temp1 	# H

cut -d ' ' -f 2 B_out_b.res >temp2 	# \Im(B^+)

# Select only the first $lines lines, corresponding to energy levels H=-1.56
# -> -1.827
# \Re(B_{in}^b), \Im(B_{in}^b)
cut -d ' ' -f 6,7 inner_ell.res | head -n $lines >temp3 

cut -d ' ' -f 3 omega_b.dat >temp4 	# \omega_{in}^b

# H, \Im(B^+), \Re(B_{in}^b), \Im(B_{in}^b), \omega_{in}^b
paste -d ' ' temp1 temp2 temp3 temp4 > B_b.dat 	
rm temp*

B_b <B_b.dat >B_b.res

# now simply plot the file B_b.res with gnuplot using B_b.plt
