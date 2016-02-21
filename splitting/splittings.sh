cut -d ' ' -f 1 ../sec1sec2/sec1sec2_notall.res  > temp1.txt  # H,
cut -d ' ' -f 4-5 ../hyper/hypers_notall.res  > temp2.txt      	# v_u[2]
cut -d ' ' -f 2 ../approxint/approxints.res  > temp3.txt   # n
cut -d ' ' -f 2-3 ../intersec/intersecs.res > temp4.txt   # p_u[2]

# Header: mu
echo "0.95387536e-3" > splittings.dat

paste -d ' ' temp1.txt temp2.txt temp3.txt temp4.txt >> splittings.dat

rm temp1.txt temp2.txt temp3.txt temp4.txt

splitting <splittings.dat >splittings.res

rm splittings.dat
