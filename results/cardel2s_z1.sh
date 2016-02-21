# Convert homoclinic point z1 (corresponding to outer splitting) from
# Cartesian to delaunay.
# Notice that z1 is on the poincare section $\Sigma_-$, so we need to use
# program cardel2s_2d.

echo "0.95387536e-3" >cardel2s_z1.dat
cut -d ' ' -f 1,5 intersec_p3.res >temp1.dat	# H, z1[0]
cut -d ' ' -f 4 ../portbp/porbits.res >temp2.dat	# z1[1]=0
paste -d ' ' temp1.dat temp2.dat >>cardel2s_z1.dat
rm temp1.dat temp2.dat
cardel2s_2d <cardel2s_z1.dat >cardel2s_z1.res
