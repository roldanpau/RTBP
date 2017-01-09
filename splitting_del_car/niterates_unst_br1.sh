# Header: mu
echo "0.95387536e-3 " > niterates_unst_br1.dat 

# H, p_u[2], z[4]
cut -d ' ' -f 1-3,9-12 ../intersec_del_car/intersecs_unst_br1.res >> niterates_unst_br1.dat 

./niterates <niterates_unst_br1.dat >niterates_unst_br1.res

#rm niterates_unst_br1.dat
