# Intersection of st mfld of p3 with horiz axis.

# This should give exactly the same result as intersec_p2.sh, i.e.
# intersection of unst mfld of p2 with horiz axis.

cut -d ' ' -f 1-3 ../sec1sec2s_p3.res > temp1.txt   # H, p[2]
cut -d ' ' -f 6-7 hyper2s_p3.res > temp2.txt      # v_s[2]
cut -d ' ' -f 3 hyper2s_p3.res > temp3.txt        # rho_s

paste -d ' ' temp1.txt temp2.txt temp3.txt > temp4.txt
tac temp4.txt >> temp5.txt

cut -d ' ' -f 2-4 approxint_p3.res > temp6.txt   # n, h_1, h_2


# Header: mu, st flag, axis line $p_x=0.0$
echo "0.95387536e-3 1 0.0" > intersec_p3.dat

paste -d ' ' temp5.txt temp6.txt >>intersec_p3.dat

rm temp*.txt
