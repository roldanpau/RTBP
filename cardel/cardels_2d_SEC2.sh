echo "0.95387536e-3" >cardels_2d_SEC2.dat
cut -d ' ' -f 1-3 ../sec1sec2/sec1sec2_2d.res >>cardels_2d_SEC2.dat
cardels_2d <cardels_2d_SEC2.dat >cardels_2d_SEC2.res
rm cardels_2d_SEC2.dat

# Convert 3-periodic point in SEC2 from Cartesian to Delaunay.