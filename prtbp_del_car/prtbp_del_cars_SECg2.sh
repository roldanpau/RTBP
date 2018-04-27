echo "0.95387536e-3 SECg2 1" >prtbp_del_cars_SECg2.dat
cut -d ' ' -f 1-4 ../cardel/cardels_SEC2.res >dels.dat
cut -d ' ' -f 1-4 ../sec1sec2/sec1sec2.res >cars.dat
paste -d ' ' dels.dat cars.dat >>prtbp_del_cars_SECg2.dat
prtbp_del_car <prtbp_del_cars_SECg2.dat >prtbp_del_cars_SECg2.res
rm dels.dat cars.dat prtbp_del_cars_SECg2.dat

# After transporting all periodic points from Cartesian to Delaunay 
# using ../cardel/cardels_SEC2.sh, we flow them to Delaunay section 
# {g=\pi} here.
