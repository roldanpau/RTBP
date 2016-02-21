# In cartesian, we obtained z1^u, a point in the unstab fundam domain of the
# periodic point $p_3$ (outer splitting).

# After transforming to delaunay, we do not obtain a point on the section
# {g=0}, and we may not obtain a point in the unst fundam domain of the
# periodic point $p_4$.

# Thus we need to iterate homoclinic point z1^u until it lies in the poincare
# section {g=0} and in the unstable fundamental domain of the point $p_4$.

echo "0.95387536e-3 1" >z1_u.dat	#mu, number of poincare iterates

cut -d ' ' -f 2-5 cardel2s_z1_u.res >>z1_u.dat	# z1_u

# WARNING! We know that we are in the unstable manifold, and we should thus
# integrate forwards using prtbpdel, NOT prtbpdel_inv.
# However, z1^u is ALMOST on the poincare section {g=0}, so we only need to
# flow it a little bit backwards to put it in the section.
# Else, we would heve to flow it a lot forwards (almost 14\pi) to have it
# again in the section, and in the unstable manifold of $p_4$.
prtbpdel_inv <z1_u.dat >z1_u.res
