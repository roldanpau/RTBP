echo "0.95387536e-3 SECg 1" >prtbpdels.dat
cut -d ' ' -f 2-5 ../cardel/cardels_2d_SEC2.res >>prtbpdels.dat
prtbpdel <prtbpdels.dat >prtbpdels.res
rm prtbpdels.dat

# After transporting all periodic points from Cartesian to Delaunay 
# using ../cardel/cardels_2d.sh, we flow them to Delaunay section 
# {g=0} here.
