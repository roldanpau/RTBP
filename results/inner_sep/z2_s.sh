# In cartesian, we obtained z2^s, a point in the stab fundam domain of the
# periodic point $p_3$ (inner splitting).

# After transforming to delaunay, we do not obtain a point on the section
# {g=0}, and we may not obtain a point in the st fundam domain of the
# periodic point $p_4$.

# Thus we need to iterate homoclinic point z2^s until it lies in the poincare
# section {g=0} and in the stable fundamental domain of the point $p_4$.

echo "0.95387536e-3 1" >z2_s.dat
cut -d ' ' -f 2-5 cardel2s_z2_s.res >>z2_s.dat	# z2_s
prtbpdel_inv <z2_s.dat >z2_s.res
