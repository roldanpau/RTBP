cut -d ' ' -f 1-4 ../sec1sec2/sec1sec2.res >cardels_SEC2.dat
cardel <cardels_SEC2.dat >cardels_SEC2.res
rm cardels_SEC2.dat

# Convert 3-periodic point in SEC2 from Cartesian to Delaunay.
