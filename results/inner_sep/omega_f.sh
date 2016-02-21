# Compute integral $\omega^f$ related to outer map.

lines=`wc -l omega_pos_f.res | cut -d ' ' -f 1`

echo $lines

# Select only the first $lines lines, corresponding to energy levels H=-1.56
# -> -1.827
cut -d ' ' -f 1 ../../portbp/porbits.res | tac | head -n $lines >temp1 	# H

head -n $lines omega_in_f.res > temp2 	# \omega_in

# H, \omega_pos, \omega_in
paste -d ' ' temp1 omega_pos_f.res temp2 > omega_f.dat 	
rm temp1 temp2 

# Compute omega_f from omega_pos and omega_in
awk -f omega_f.awk <omega_f.dat >omega_f.res 

# now simply plot the file omega_f.res with gnuplot using omega_f.plt
