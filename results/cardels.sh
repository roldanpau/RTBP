# Convert periodic point p3 (located at apohelion) from Cartesian to delaunay

echo "0.95387536e-3" >cardels.dat
cut -d ' ' -f 1,3,4 ../portbp/porbits.res >>cardels.dat
cardels_2d <cardels.dat >cardels.res
