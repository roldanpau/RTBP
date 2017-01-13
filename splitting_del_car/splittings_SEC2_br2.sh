datfile=splittings_SEC2_br2.dat
resfile=splittings_SEC2_br2.res

cut -d ' ' -f 1 ../sec1sec2/sec1sec2.res  > temp1.txt   # H,
cut -d ' ' -f 2 niterates_SEC2_br2.res > temp2.txt      # sec
cut -d ' ' -f 4-5 ../hyper/hypers_SEC1.res  > temp3.txt # v_u[2] in SEC1
cut -d ' ' -f 4-5 ../hyper/hypers.res  > temp4.txt      # v_u[2] in SEC2
cut -d ' ' -f 3 niterates_SEC2_br2.res > temp5.txt      # n
cut -d ' ' -f 2-3 ../intersec_del_car/intersecs_unst_SEC2_br2.res > temp6.txt   # p_u[2]

# Header: mu
echo "0.95387536e-3 " > $datfile

paste -d ' ' temp1.txt temp2.txt temp3.txt temp4.txt temp5.txt temp6.txt >> $datfile

rm temp*.txt

./splitting_del_car < $datfile > $resfile

rm $datfile
