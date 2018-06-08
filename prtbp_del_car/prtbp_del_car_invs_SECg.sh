DATFILE=prtbp_del_car_invs_SECg.dat
RESFILE=prtbp_del_car_invs_SECg.res

echo "0.95387536e-3 SECg 1" >$DATFILE
cut -d ' ' -f 1-4 ../cardel/cardels_SEC2.res >dels.dat
cut -d ' ' -f 1-4 ../sec1sec2/sec1sec2.res >cars.dat
paste -d ' ' dels.dat cars.dat >>$DATFILE
prtbp_del_car_inv <$DATFILE >$RESFILE
rm dels.dat cars.dat $DATFILE

# After transporting all periodic points from Cartesian to Delaunay 
# using ../cardel/cardels_SEC2.sh, we flow them to Delaunay section 
# {g=0} here.
# We flow them using prtbp_del_car_inv just to test this program.
