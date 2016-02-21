tail -n +2 approxint_p3.dat | cut -d ' ' -f 1-6 - > temp1.txt	# H, p[2], v_u[2], rho_u
cut -d ' ' -f 2-4 approxint_p3_m2h.res > temp2.txt		# n, h_1, h_2

# Header: mu, unst flag, axis line $p_x=-2h$
echo "0.95387536e-3 0 -2.e-4" > intersec_p3_m2h.dat

paste -d ' ' temp1.txt temp2.txt >>intersec_p3_m2h.dat

rm temp*.txt

intersec <intersec_p3_m2h.dat >intersec_p3_m2h.res
