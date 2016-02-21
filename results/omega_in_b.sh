# Compute partial shift $\omega_{in}^b$. This is needed to obtain the outer
# map.

# We output lines with INCREASING energy H=-1.56 -> -2.04

echo "0.95387536e-3 0" >omega_in_b.dat	# mu, bInner flag

# reverse lines, to get INCREASING energy values
cut -d ' ' -f 1-4 prtbpdel_p3.res | tac >> omega_in_b.dat # periodic point p3
nohup inner_circ < omega_in_b.dat > omega_in_b.res 2> omega_in_b.err &
