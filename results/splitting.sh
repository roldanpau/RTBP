cut -d ' ' -f 1 sec1sec2s_p3.res > temp1.txt   # H
cut -d ' ' -f 4-5 hyper2s_p3.res > temp2.txt   # v_u[2]

paste -d ' ' temp1.txt temp2.txt > temp3.txt
tac temp3.txt >> temp4.txt

cut -d ' ' -f 2 approxint_p3.res > temp5.txt   # n
cut -d ' ' -f 2,3 intersec_p3.res > temp6.txt   # p_u[2]

# Header: mu
echo "0.95387536e-3" > splitting.dat

paste -d ' ' temp4.txt temp5.txt temp6.txt >>splitting.dat

rm temp*.txt
