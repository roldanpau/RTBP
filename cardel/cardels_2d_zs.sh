echo "0.95387536e-3" >cardels_2d_zs.dat

# Note: zs[2] = -zu[2]
cut -d ' ' -f 1-3 ../intersec/intersecs.res | awk '{print $1, $2, -$3}' >>cardels_2d_zs.dat	# H, zs[2]

cardels_2d <cardels_2d_zs.dat >cardels_2d_zs.res
rm cardels_2d_zs.dat
