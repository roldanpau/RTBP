# Compute function $B_{out}^{b,+}$ associated to outer map of elliptic
# problem.
#
# Actually, we only compute the integral called $B^+$ in my notes. This is
# enough, because 
#    B_{out}^{b,+} = -\mu B^+ + \mu C^+$,
# where
#    re(C^+)=re(B^+), and im(C^+)=-im(B^+).

lines=`wc -l omega_pos_b.res | cut -d ' ' -f 1`

echo $lines

# Select only the first $lines lines, corresponding to energy levels H=-1.56
# -> -1.827
cut -d ' ' -f 1 ../portbp/porbits.res | tac | head -n $lines >temp1 	# H=-1.56 -> -2.04

# Select only the first $lines lines, corresponding to energy levels H=-1.56
# -> -1.827
cut -d ' ' -f 1-4 prtbpdel_p4.res | tac | head -n $lines >temp2 	# p4

cut -d ' ' -f 1-4 z1_u.res >temp3 	# z1_u
cut -d ' ' -f 2 approxint_p3.res >temp4 # N

echo "0.95387536e-3" > B_out_b.dat

# H, p4, z1_u, \omega_pos, N
paste -d ' ' temp1 temp2 temp3 omega_pos_b.res temp4 >> B_out_b.dat
rm temp*

nohup outer_ell <B_out_b.dat >B_out_b.res 2>B_out_b.err &
