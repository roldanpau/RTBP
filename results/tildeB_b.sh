# Compute function $\tilde B^{b,+)$ related to averaged outer map.

lines=`wc -l omega_b.dat | cut -d ' ' -f 1`

echo $lines

# Select only the first $lines lines, corresponding to energy levels H=-1.56
# -> -1.827
# H, T_H
cut -d ' ' -f 1,2 ../portbp/porbits.res | tac | head -n $lines >temp1 	

cut -d ' ' -f 2 omega_b.res >temp2	# omega^b

# re(A_1^+), im(A_1^+)
cut -d ' ' -f 2,3 inner_ell.res | head -n $lines >temp3

# re(B^{b,+}), im(B^{b,+})
cut -d ' ' -f 2,3 B_b.res >temp4

# H, T_H, omega^b, re(A_1^+), im(A_1^+), re(B^{b,+}), im(B^{b,+})
paste -d ' ' temp1 temp2 temp3 temp4 > tildeB_b.dat 	
rm temp*

# Finally, compute \tilde B^b using program tOmega in directory tOmega
