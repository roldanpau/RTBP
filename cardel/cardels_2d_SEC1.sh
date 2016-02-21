echo "0.95387536e-3" >cardels_2d_SEC1.dat
cut -d ' ' -f 1,3,5 ../portbp/porbits.res >>cardels_2d_SEC1.dat
cardels_2d <cardels_2d_SEC1.dat >cardels_2d_SEC1.res
rm cardels_2d_SEC1.dat

# Convert 3-periodic point p_0 in SEC1 from Cartesian to Delaunay.
# Since p_0 is in the symmetry line, in Delaunay it should have $\ell=0$.
