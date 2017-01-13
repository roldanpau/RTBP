datfile=niterates_SEC1_br2.dat
resfile=niterates_SEC1_br2.res

# Header: mu
echo "0.95387536e-3 " > $datfile

# H, p_u[2], z[4]
cut -d ' ' -f 1-3,9-12 ../intersec_del_car/intersecs_unst_SEC1_br2.res >> $datfile

./niterates <$datfile >$resfile

rm $datfile
