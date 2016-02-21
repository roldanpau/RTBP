cut -d ' ' -f 1,5 intersec_p3_p1h.res > temp1.txt   # H, x_u^1
cut -d ' ' -f 5 intersec_p2_p1h.res > temp2.txt   	# x_s^1
cut -d ' ' -f 5 intersec_p3_m1h.res > temp3.txt # x_u^{-1}
cut -d ' ' -f 5 intersec_p2_m1h.res > temp4.txt # x_s^{-1}
cut -d ' ' -f 5 intersec_p3_p2h.res > temp5.txt  # x_u^2
cut -d ' ' -f 5 intersec_p2_p2h.res > temp6.txt  # x_s^2
cut -d ' ' -f 5 intersec_p3_m2h.res > temp7.txt # x_u^{-2}
cut -d ' ' -f 5 intersec_p2_m2h.res > temp8.txt # x_s^{-2}

# Header: mu
echo "0.95387536e-3" > axis.dat

paste -d ' ' temp1.txt temp2.txt temp3.txt temp4.txt temp5.txt temp6.txt temp7.txt temp8.txt >>axis.dat
#paste -d ' ' temp1.txt temp2.txt temp3.txt temp4.txt >>axis.dat

rm temp*.txt

./axis <axis.dat >axis.res
