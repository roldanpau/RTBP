# Convert homoclinic point z2 (corresponding to inner splitting) from
# Cartesian to delaunay.
# Notice that z2 is on the poincare section $\Sigma_-$, so we need to use
# program cardel2s_2d.

echo "0.95387536e-3" >cardel2s_z2.dat
cut -d ' ' -f 1,5 intersec_p2.res >temp1.dat	# H, z2[0]
cut -d ' ' -f 4 ../../portbp/porbits.res >temp2.dat	# z2[1]=0
paste -d ' ' temp1.dat temp2.dat >>cardel2s_z2.dat
rm temp1.dat temp2.dat

# Remember to delete last lines of cardel2s_z2.dat, and then do:
#cardel2s_2d <cardel2s_z2.dat >cardel2s_z2.res

#Finally, remember to zero g component in cardel2s_z2.res (if necessary).
