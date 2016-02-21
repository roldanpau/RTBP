# Test computed value of \omega_-^f by checking that homoclinic orbit
# \gamma^f converges to periodic orbit \gamma^4

lines=`wc -l omega_pos_f.res | cut -d ' ' -f 1`

echo $lines

cut -d ' ' -f 1 ../../portbp/porbits.res | tac >temp1 	# H=-1.56 -> -2.04
# Select only the first $lines lines, corresponding to energy levels H=-1.56
# -> -1.827
head -n $lines temp1 > temp2 		# H

cut -d ' ' -f 1-4 ../prtbpdel_p4.res | tac >temp3
# Select only the first $lines lines, corresponding to energy levels H=-1.56
# -> -1.827
head -n $lines temp3 > temp4 		# p3

cut -d ' ' -f 1-4 z2_s.res >temp5 	# z2_s
cut -d ' ' -f 2 approxint_p3.res >temp6 # N

echo "0.95387536e-3" > outer_circ_test_neg.dat

# Since omega_neg = -omega_pos, pass omega_pos, and program
# outer_circ_test_neg will negate this value.

# H, p4, z2_s, N, \omega_pos
paste -d ' ' temp2 temp4 temp5 temp6 omega_pos_f.res >> outer_circ_test_neg.dat
rm temp*

outer_circ_test_neg <outer_circ_test_neg.dat >outer_circ_test_neg.res
