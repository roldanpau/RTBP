echo "0.95387536e-3" >cardels_2d_zf.dat
cut -d ' ' -f 1,5 ../intersec/intersecs.res >temp1.dat
awk '{print $1, $2, 0}' temp1.dat >>cardels_2d_zf.dat
cardels_2d <cardels_2d_zf.dat >cardels_2d_zf.res
rm temp1.dat cardels_2d_zf.dat
