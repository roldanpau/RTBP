tail -n +2 approxint_p2.dat | cut -d ' ' -f 1-6 - > temp1.txt	# H, p[2], v_s[2], rho_s
cut -d ' ' -f 2-4 approxint_p2_m2h.res > temp2.txt		# n, h_1, h_2

# Header: mu, st flag, axis line $p_x=-2h$
echo "0.95387536e-3 1 -2.e-4" > intersec_p2_m2h.dat

paste -d ' ' temp1.txt temp2.txt >>intersec_p2_m2h.dat

rm temp*.txt

intersec <intersec_p2_m2h.dat >intersec_p2_m2h.res
