# Convert homoclinic point z2^s in stable fundamental domain (corresponding
# to inner splitting) from Cartesian to delaunay.
# Notice that z2^s is on the poincare section $\Sigma_-$, so we need to use
# program cardel2s_2d.

# Warning! the output of this program is a point in delaunay, but it is not
# necessarily on the poincare section {g=0}! Thus we need to further process
# it using z2_s.sh

echo "0.95387536e-3" >cardel2s_z2_s.dat
cut -d ' ' -f 1,2-3 intersec_p3.res >>cardel2s_z2_s.dat	# H, z2_s
cardel2s_2d <cardel2s_z2_s.dat >cardel2s_z2_s.res
