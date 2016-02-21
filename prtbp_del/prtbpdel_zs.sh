echo "0.95387536e-3 SEC1 1" >prtbpdel_zs.dat
cut -d ' ' -f 2-5 ../cardel/cardels_2d_zs.res >>prtbpdel_zs.dat
prtbpdel_inv <prtbpdel_zs.dat >prtbpdel_zs.res
rm prtbpdel_zs.dat
