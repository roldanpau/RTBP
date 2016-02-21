echo "0.95387536e-3" >cardels_2d_zu.dat
cut -d ' ' -f 1-3 ../intersec/intersecs.res >>cardels_2d_zu.dat	# H, zu[2]
cardels_2d <cardels_2d_zu.dat >cardels_2d_zu.res
rm cardels_2d_zu.dat
