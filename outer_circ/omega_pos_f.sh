# Compute integral $\omega_pos^f$ related to outer map.

echo "0.95387536e-3" >omega_pos_f.dat   # mu

cut -d ' ' -f 2 ../portbp/porbits_notall.res > temp1    # period T

# z_s: preimage of homoclinic point z
cut -d ' ' -f 1-4 ../prtbp_del/prtbpdel_zs.res > temp2

# M: number of iterations to reach homoclinic point z from z_s
cut -d ' ' -f 2 ../approxint/approxints.res > temp3

paste -d ' ' temp1 temp2 temp3 >>omega_pos_f.dat
rm temp1 temp2 temp3

nohup outer_circ < omega_pos_f.dat > omega_pos_f.res 2> omega_pos_f.err &
rm omega_pos_f.dat
