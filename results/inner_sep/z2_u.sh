# In cartesian, we obtained z2^u, a point in the unstab fundam domain of the
# periodic point $p_2$ (inner splitting).

# After transforming to delaunay, we do not obtain a point on the section
# {g=0}, and we may not obtain a point in the unst fundam domain of the
# periodic point $p_3$.

# Thus we need to iterate homoclinic point z2^u until it lies in the poincare
# section {g=0} and in the unstable fundamental domain of the point $p_3$.

echo "0.95387536e-3 1" >z2_u.dat
cut -d ' ' -f 2-5 cardel2s_z2_u.res >>z2_u.dat	# z2_u
prtbpdel <z2_u.dat >z2_u.res
