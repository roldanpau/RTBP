cut -d ' ' -f 1-3 ../sec1sec2s_p2.res > temp1.txt   # H, p[2]
cut -d ' ' -f 4-5 hyper2s_p2.res > temp2.txt      # v_u[2]
cut -d ' ' -f 2 hyper2s_p2.res > temp3.txt        # rho_u

# Header: mu, k, unstable flag, axis line $p_x=a$
echo "0.95387536e-3 6 0 0.0" > approxint_p2.dat

paste -d ' ' temp1.txt temp2.txt temp3.txt >> approxint_p2.dat

rm temp1.txt temp2.txt temp3.txt
