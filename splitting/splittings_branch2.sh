cut -d ' ' -f 1 ../sec1sec2/sec1sec2_notall.res  > temp1.txt  # H,
cut -d ' ' -f 6-7 ../hyper/hypers_notall.res  > temp2.txt      	# v_s[2]
cut -d ' ' -f 2 ../approxint/approxints_branch2.res  > temp3.txt   # n
cut -d ' ' -f 2-3 ../intersec/intersecs_branch2.res > temp4.txt   # p_s[2]

# Header: mu
echo "0.95387536e-3" > splittings_branch2.dat

paste -d ' ' temp1.txt temp2.txt temp3.txt temp4.txt >> splittings_branch2.dat

rm temp1.txt temp2.txt temp3.txt temp4.txt

splitting <splittings_branch2.dat >splittings_branch2.res

rm splittings_branch2.dat
