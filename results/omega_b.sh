# Compute integral $\omega^b$ related to outer map.

lines=`wc -l omega_pos_b.res | cut -d ' ' -f 1`

echo $lines

# Select only the first $lines lines, corresponding to energy levels H=-1.56
# -> -1.827
cut -d ' ' -f 1 ../portbp/porbits.res | tac | head -n $lines >temp1 	# H

head -n $lines omega_in_b.res > temp2 	# \omega_in

# H, \omega_pos, \omega_in
paste -d ' ' temp1 omega_pos_b.res temp2 > omega_b.dat 	
rm temp1 temp2 

awk -f omega_b.awk <omega_b.dat >omega_b.res

# now simply plot the file omega_b.res with gnuplot using omega_fb.plt
