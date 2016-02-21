# Compute integral $\omega_pos^b$ related to outer map.

lines=`wc -l z1_u.res | cut -d ' ' -f 1`

echo $lines

# We output lines with INCREASING energy H=-1.56 -> -2.04

echo "0.95387536e-3" >omega_pos_b.dat	# mu

# reverse lines, to get INCREASING energy values
# Select only the first $lines lines, corresponding to energy levels H=-1.56
# -> -1.827
cut -d ' ' -f 2 ../portbp/porbits.res | tac | head -n $lines > temp1	# period T

# z1_u.res is already in INCREASING energy values
# z_u: preimage of homoclinic point z
cut -d ' ' -f 1-4 z1_u.res > temp2

# approxint_p3.res is already in INCREASING energy values
# M: number of iterations to reach homoclinic point z from z_u
cut -d ' ' -f 2 approxint_p3.res > temp3

paste -d ' ' temp1 temp2 temp3 >>omega_pos_b.dat
rm temp1 temp2 temp3

nohup outer_circ < omega_pos_b.dat > omega_pos_b.res 2> omega_pos_b.err &
