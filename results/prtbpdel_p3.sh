# Take periodic point (located at apohelion) to Poincare section {g=0} by
# flowing backwards. This gives periodic point p3 in Delaunay.

echo "0.95387536e-3 1" >prtbpdel_p3.dat
cut -d ' ' -f 2-5 cardels.res >>prtbpdel_p3.dat
prtbpdel_inv <prtbpdel_p3.dat >prtbpdel_p3.res
