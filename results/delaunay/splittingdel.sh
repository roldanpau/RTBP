cut -d ' ' -f 1,2 sec1sec2s_p0.res > temp1.txt   # e,H
cut -d ' ' -f 5-6 hyperdels_p0.res > temp2.txt   # v_u[2]

cut -d ' ' -f 3 approxintdel_p0.res > temp3.txt   # n
cut -d ' ' -f 3,4 intersecdel_p0.res > temp4.txt   # p_u[2]

# Header: mu
echo "0.95387536e-3" > splittingdel.dat

paste -d ' ' temp1.txt temp2.txt temp3.txt temp4.txt >>splittingdel.dat

rm temp*.txt
