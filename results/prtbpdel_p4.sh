# Take periodic point (located at apohelion) to Poincare section {g=0} by
# flowing forwards. This gives periodic point p4 in Delaunay.

echo "0.95387536e-3 1" >prtbpdel_p4.dat
cut -d ' ' -f 2-5 cardels.res >>prtbpdel_p4.dat
prtbpdel <prtbpdel_p4.dat >prtbpdel_p4.res
