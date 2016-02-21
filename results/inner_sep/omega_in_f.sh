# Compute partial shift $\omega_{in}^f$. This is needed to obtain the outer
# map.

# We output lines with INCREASING energy H=-1.56 -> -2.04

echo "0.95387536e-3 1" >omega_in_f.dat	# mu, bInner flag

# reverse lines, to get INCREASING energy values
cut -d ' ' -f 1-4 ../prtbpdel_p4.res | tac >> omega_in_f.dat # periodic point p4
nohup inner_circ < omega_in_f.dat > omega_in_f.res 2> omega_in_f.err &
