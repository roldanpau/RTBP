# Compute function $\tilde B^{f,+)$ related to averaged outer map.

lines=`wc -l omega_pos_f.res | cut -d ' ' -f 1`

echo $lines

# Select only the first $lines lines, corresponding to energy levels H=-1.56
# -> -1.827
# H, T_H
cut -d ' ' -f 1,2 ../../portbp/porbits.res | tac | head -n $lines >temp1 	

cut -d ' ' -f 2 omega_f.res >temp2	# omega^f

# re(A_1^+), im(A_1^+)
cut -d ' ' -f 2,3 ../inner_ell.res | head -n $lines >temp3

# re(B^{f,+}), im(B^{f,+})
cut -d ' ' -f 2,3 B_f.res >temp4

# H, T_H, omega^f, re(A_1^+), im(A_1^+), re(B^{f,+}), im(B^{f,+})
paste -d ' ' temp1 temp2 temp3 temp4 > tildeB_f.dat 	
rm temp*

# Finally, compute \tilde B^f with program tOmega in directory tOmega.
