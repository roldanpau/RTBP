# Compute integral $\omega^f$ related to outer map.

cut -d ' ' -f 1 ../portbp/porbits_notall.res >temp1 	# H

# H, \omega_neg, \omega_in
paste -d ' ' temp1 omega_neg_f.res ../inner_circ/omega_in_f.res > omega_f.dat 	
rm temp1

#awk -f omega_f.awk <omega_f.dat >omega_f.res
#rm omega_f.dat

# now simply plot the file omega_f.res with gnuplot using omega_fb.plt
