cut -d ' ' -f 1-3 sec1sec2s_p2.res > temp1.txt   # H, p[2]
cut -d ' ' -f 6-7 hyper2s_p2.res > temp2.txt      # v[2]
cut -d ' ' -f 3 hyper2s_p2.res > temp3.txt        # rho

# Header: mu, k, stable flag, axis line $p_x=a$
echo "0.95387536e-3 6 1 0.0" > approxint_p2.dat

paste -d ' ' temp1.txt temp2.txt temp3.txt >> approxint_p2.dat

rm temp1.txt temp2.txt temp3.txt
